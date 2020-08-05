"""
This file contains a  class for extracting
intersections and differences of networks.
These are returned as graphml files. """

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import sys
import os
import networkx as nx
import logging.handlers
from masq.scripts.utils import ParentConnection

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


def start_metastats(level, networks=None,
                    weight=True,
                    config='database.ini',
                    host=None, database=None,
                    username=None, password=None):
    """
    Can read a single network or all networks in a folder.
    These are then imported into the sqlite database.
    Unless a mapping file is supplied,
    the filename is used to generate the network name.
    The network name should match the BIOM file,
    otherwise metadata (e.g. node taxonomy) cannot be queried.

    :param set: Type of set to extract
    :param networks: Networks to extract set from
    :param config: Location of file with database parameters.
    :param host: Database address.
    :param database: Name of PostgreSQL database.
    :param username: Username for PostgreSQL database.
    :param password: Password of PostgreSQL database.
    :return:
    """
    conn = MetaConnection(config, host, database, username, password)
    if not networks:
        networks = conn.get_networks()
    tax_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
    level_id = tax_list.index(level.capitalize())
    for level in range(0, level_id + 1):
        logger.info("Checking " + tax_list[level] + " level...")
        networks = MetaConnection.agglomerate_networks(level=tax_list[level],
                                                       weight=weight, networks=networks)


class MetaConnection(ParentConnection):
    """
    Initializes a connection to the PostgreSQL database.
    This connection object contains methods for aggregating networks
    by taxonomic level.
    """
    # inherits init from parent
    def agglomerate_networks(self, level=None, weight=True, networks=None):
        """
        Agglomerates to specified taxonomic level, or, if no level is specified,
        over all levels. Edges are agglomerated based on similarity
        at the specified taxonomic level. If 'weight' is set to True,
        edges are only agglomerated if their weight matches.
        The stop condition is the length of the pair list;
        as soon as no pair meets the qualification, agglomeration is terminated.
        By default, agglomeration is done separately
        per network in the database, so each network gets an agglomerated version.

        The networks parameter can be both a dict and a list.
        If it is a dict, the keys are the new network names, the values the old names.

        Pseudocode representation:
        1. Duplicate networks
        2. For each edge pair (taxon-level)-taxon-taxon-(taxon-level)
            3. Create new edge
            4. Delete edge pair

        :param networks: List or dict of networks to extract agglomerated network from.
        :param level: Taxonomic level matching taxonomic assignments in the PostgreSQL database
        :param weight: If true, an edge is counted separately if the sign of the weight is different
        :return:
        """
        tax_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
        new_networks = dict()
        # we create a copy of the original network
        for network in networks:
            previous_network = network
            if type(networks) == list:
                new_name = level + '_' + network
            else:
                new_name = level + '_' + '_'.join(network.split('_')[1:])
            new_networks[new_name] = network
            # first check if lower-level network exists
            # if there were no pairs, it might not have been copied
            network_query = "SELECT networkID from networks WHERE networkID = %s"
            hit = self.value_query(query=network_query,
                                   values=(network,),
                                   fetch=True)
            i = tax_list.index(level)
            while len(hit) == 0:
                previous_network = '_'.join(network.split('_')[1:])
                if i > 0:
                    previous_network = tax_list[i] + '_' + previous_network
                    i -= 1
                hit = self.value_query(query=network_query,
                                       values=(previous_network,),
                                       fetch=True)
            # if there are no pairs at all, no need to copy network
            # possible with nodes that do not have large taxonomy
            testpair = self.get_pairlist(level, weight, previous_network)
            if len(testpair) == 0:
                testpair = self.get_taxlist(level, previous_network)
            if not len(testpair) == 0:
                logger.info("Copying " + previous_network + "...")
                self.copy_network(previous_network, new_name)
            else:
                new_networks[new_name] = None
        try:
            for network in new_networks:
                if new_networks[network]:
                    logger.info("Agglomerating " + network + "...")
                    stop_condition = False
                    while not stop_condition:
                        # limit is necessary to prevent excessively long queries
                        pairs = self.get_pairlist(level=level, weight=weight, network=network)
                        if pairs:
                            pair = pairs[0]['p']
                            self.agglomerate_pair(pair, level=level, weight=weight, network=network)
                        else:
                            stop_condition = True
                    stop_condition = False
                    while not stop_condition:
                        # after agglomerating edges
                        # taxa with same taxonomic assignments should be merged
                        # this rewires the network
                        tax_nodes = self.get_taxlist(level=level, network=network)
                        if tax_nodes:
                            tax_nodes = tax_nodes[0]['p']
                            self.agglomerate_taxa(tax_nodes, level=level, weight=weight)
                        else:
                            stop_condition = True
                    num = self.query("MATCH (n:Network {name: '" + network + "'})-"
                                     "[r:PART_OF]-() RETURN count(r) as count")
                    logger.info("The agglomerated network " + network +
                                " contains " + str(num) + " edges.")
        except Exception:
            logger.error("Could not agglomerate edges to higher taxonomic levels. \n", exc_info=True)
        return new_networks

    def get_pairlist(self, level, weight, network):
        """
        Returns a single pair of edges.
        A pair is defined as two edges linked to taxonomic nodes
        that have the same taxonomic assignment at the specified level,
        e.g. Nitrobacter-edge-Nitrosomonas.

        Pseudo-code for SQL query:
        For each row of the edge table
        Join source - taxonomy table
        Join target - taxonomy table
        Group by edges with same source taxonomy and target taxonomy
        Return source, target IDs

        Problem: source and target cannot be swapped internally
        for the GROUP BY statement. At least not as fas as I know.
        Solution: do query on copy of edges, with each edge copied in reverse.
        The taxonomy join should then show edges added with reverse taxonomy.

        :param level: Taxonomic level to identify a pair.
        :param weight: if True, specifies that edge weights should have the same sign.
        :param network: Name of network that the pairs should belong to
        :return: List containing results of Neo4j transaction
        """
        # first, modify the table so that
        # edge taxonomy at the level of interest is in alphabetical order
        copy_query = "CREATE TABLE copy AS " \
                     "SELECT source, target, weight FROM edges " \
                     "WHERE edges.networkID = %s;"
        self.value_query(copy_query, (network,))
        copy_query = "INSERT INTO copy (source, target, weight) " \
                     "SELECT target, source, weight FROM copy;"
        self.value_query(copy_query, (network,))
        agglom_query = "SELECT string_agg(source::varchar, ',') as source," \
                       "string_agg(target::varchar, ',') as target, p." + level + \
                       " as source, q." + level + " as target FROM copy as e " \
                       "JOIN taxonomy as p ON e.source = p.taxon " \
                       "JOIN taxonomy as q on e.target = q.taxon " \
                       "GROUP BY p." + level + ", q." + level + \
                       " HAVING COUNT(*) > 1 LIMIT 1;"
        results = self.value_query(agglom_query, values=(network,), fetch=True)
        self.query("DROP TABLE copy;")
        sources = results[0][0].split(',')
        targets = results[0][1].split(',')
        # need to check if sources and targets need to be swapped
        check_query = "SELECT source, target from edges WHERE networkID = %s " \
                      "AND source=%s AND target=%s"
        for i in range(len(sources)):
            checks = self.value_query(check_query, values=(network, sources[i], targets[i]), fetch=True)
            if len(checks) == 0:
                sources[i] = results[0][1].split(',')[i]
                targets[i] = results[0][0].split(',')[i]
        return sources, targets
