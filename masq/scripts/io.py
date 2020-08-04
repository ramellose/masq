"""
This file contains a  class for writing a network file
to a PostgreSQL database.
 Both a summary network table are edited,
 as well as an edge list table. """

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


def import_networks(location, mapping=None, sources=None,
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

    :param location: Location of network file(s)
    :param mapping: Dictionary with BIOM filenames as keys, new names as values
    :param sources: Dictionary with network filenames as keys, BIOM names as values
    :param config: Location of file with database parameters.
    :param host: Database address.
    :param database: Name of PostgreSQL database.
    :param username: Username for PostgreSQL database.
    :param password: Password of PostgreSQL database.
    :return:
    """
    conn = IoConnection(config, host, database, username, password)
    if os.path.isdir(location):
        for y in os.listdir(location):
            network = _read_network_extension(location + '/' + y)
            name = y.split(".")[0]
            source = name
            if mapping:
                name = mapping[name]
                source = name
            if source:
                source = sources[name]
            conn.add_network(network, name=name, study=source)
    else:
        network = _read_network_extension(location)
        name = location.split('/')[-1]
        name = name.split('\\')[-1]
        name = name.split(".")[0]
        source = name
        if mapping:
            name = mapping[name]
            source = name
        if sources:
            source = sources[name]
        conn.add_network(network, name=name, study=source)


class IoConnection(ParentConnection):
    """
    Initializes a connection to the sqlite3 database.
    This connection object contains methods for converting NetworkX
    objects to rows in the summary network table
    and rows in the edges table.
    """
    # inherits init from parent
    def add_network(self, network, name, study):
        """
        Takes a networkx object and writes this to the sqlite3 database.
        :param network: NetworkX object
        :param name: Network name
        :param study: Study ID (needs to match a biom ID)
        :return:
        """
        node_num = len(network.nodes)
        edge_num = len(network.edges)
        network_values = name, study, node_num, edge_num
        self.add_network_node(network_values)
        edge_values = list()
        for edge in network.edges:
            if 'weight' in network.edges[edge]:
                edge_values.append((name, edge[0], edge[1], network.edges[edge]['weight']))
            else:
                edge_values.append((name, edge[0], edge[1], None))
        self.add_edge(edge_values)
        logger.info("Uploaded network data for " + name + ".\n")

    def add_network_node(self, values):
        """
        Adds rows to the networks table in the PostgreSQL database.

        :param values: List of tuples for bioms table
        :return:
        """
        network_query = "INSERT INTO networks(networkID, studyID,node_num,edge_num) " \
                        "VALUES (%s, %s,%s,%s)"
        self.value_query(network_query, values)

    def add_edge(self, values):
        """
        Adds rows to the edges table in the PostgreSQL database.

        :param values: List of tuples for edges table
        :return:
        """
        edge_query = "INSERT INTO edges (networkID,source,target,weight) " \
                     "VALUES (%s,%s,%s,%s)"
        self.value_query(edge_query, values)

    def export_network(self, name):
        """
        Extracts all edges belonging to a specific network,
        and returns these as a networkx object.
        :param name: Network name
        :return:
        """
        network_query = "SELECT * from edges WHERE " \
                        "networkID='" + name + "';"
        edges = self.query(network_query, fetch=True)
        network = nx.Graph()
        for edge in edges:
            network.add_edge(edge[1], edge[2], weight=edge[3])
        return network


def _read_network_extension(filename):
    """
    Given a filename with a specific extension,
    this function calls the correct function to read the file.

    :param filename: Complete filename.
    :return: NetworkX object
    """
    extension = filename.split(sep=".")
    extension = extension[len(extension) - 1]
    network = None
    if extension == 'graphml':
        network = nx.read_graphml(filename)
    elif extension == 'txt':
        network = nx.read_weighted_edgelist(filename)
    elif extension == 'gml':
        network = nx.read_gml(filename)
    else:
        logger.warning('Ignoring file with wrong format.', exc_info=True)
        network = False
    if network:
        try:
            if 'name' in network.nodes[list(network.nodes)[0]]:
                if network.nodes[list(network.nodes)[0]]['name'] != list(network.nodes)[0]:
                    network = nx.relabel_nodes(network, nx.get_node_attributes(network, 'name'))
        except IndexError:
            logger.warning('One of the imported networks contains no nodes.', exc_info=True)
    return network