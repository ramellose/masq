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


def extract_sets(path, set, networks=None,
                 size=None, weight=True,
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
    conn = SetConnection(config, host, database, username, password)
    if not networks:
        networks = conn.get_networks()
    if set == 'intersection':
        if not size:
            size = 1
        elif size*len(networks) > len(networks):
            logger.error("Cannot extract intersection for more networks "
                         "than the number of networks in the database. ")
        g = conn.get_intersection(networks, number=int((len(networks)*size)),
                                  weight=weight)
    elif set == 'difference':
        g = conn.get_difference(networks, weight=weight)
    elif set == 'union':
        g = conn.get_union(networks)
    nx.write_graphml(g, path + "//" + set + ".graphml")


class SetConnection(ParentConnection):
    """
    Initializes a connection to the PostgreSQL database.
    This connection object contains methods for converting NetworkX
    objects to rows in the summary network table
    and rows in the edges table.
    """
    # inherits init from parent
    def get_intersection(self, networks, number, weight=True):
        """
        :param networks: List of networks to extract intersection from
        :param number: Number of networks where an edge has to occur
        :param weight: If true, an edge is counted separately if the sign of the weight is different
        :return:
        """
        # first, extract all edges that occur more than once
        # these can have edges with different weight
        if weight:
            set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                        "source, target, SIGN(weight), COUNT(*) " \
                        "FROM edges " \
                        "WHERE networkID IN " + str(tuple(networks)) + \
                        " GROUP BY source, target, SIGN(weight) " \
                        "HAVING COUNT(*) > " + str(number-1)
        else:
            set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                        "source, target, string_agg(weight::varchar, ',') as weights, COUNT(*) " \
                        "FROM edges " \
                        "WHERE networkID IN " + str(tuple(networks)) + \
                        " GROUP BY source, target " \
                        "HAVING COUNT(*) > " + str(number-1)
        set_result = self.query(set_query, fetch=True)
        g = _convert_network(set_result)
        logger.info("Extracted intersection across " + str(number) + " networks...\n")
        return g

    def get_difference(self, networks, weight=True):
        """
        :param networks: List of networks to extract intersection from
        :param weight: If true, an edge is counted separately if the sign of the weight is different
        :return:
        """
        if weight:
            set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                        "source, target, SIGN(weight), COUNT(*) " \
                        "FROM edges " \
                        "WHERE networkID IN " + str(tuple(networks)) + \
                        " GROUP BY source, target, SIGN(weight) " \
                        "HAVING COUNT(*) = 1 "
        else:
            set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                        "source, target, string_agg(weight::varchar, ',') as weights, COUNT(*) " \
                        "FROM edges " \
                        "WHERE networkID IN " + str(tuple(networks)) + \
                        " GROUP BY source, target " \
                        "HAVING COUNT(*) = 1 "
        set_result = self.query(set_query, fetch=True)
        g = _convert_network(set_result)
        logger.info("Extracted difference...\n")
        return g

    def get_union(self, networks):
        """
        :param network: NetworkX object
        :param name: Network name
        :param study: Study ID (needs to match a biom ID)
        :return:
        """
        set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                    "source, target, string_agg(weight::varchar, ',') as weights, COUNT(*) " \
                    "FROM edges " \
                    "WHERE networkID IN " + str(tuple(networks)) + \
                    " GROUP BY source, target;"
        set_result = self.query(set_query, fetch=True)
        g = _convert_network(set_result)
        logger.info("Extracted union...\n")
        return g

    def get_networks(self):
        """
        Gets the network names from the database,
        in case these are not given by the user.
        :return: List with network names
        """
        network_query = "SELECT networkID from networks;"
        networks = self.query(network_query, fetch=True)
        networks = [x[0] for x in networks]
        return networks


def _convert_network(edge_list):
    """
    Takes a SQL output with edges (and weights).
    With the edge list having a length of 5,
    the 4th value is assumed to be the edge sign / weight.
    If the query was carried out without grouping by weight,
    weight is a string variable generated from aggregated weights.

    WARNING: networkx does not support edges with multiple weights.
    So the edge may be overwritten.
    :param edge_list:
    :return: Networkx graph
    """
    g = nx.Graph()
    if len(edge_list) > 0:
        for edge in edge_list:
            g.add_edge(edge[1], edge[2], source=edge[0], weight=edge[3])
    return g