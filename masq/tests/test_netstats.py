"""
This file contains functions for testing functions
and the ParentConnection class in the utils.py file.

The file first connects to a simple PostgreSQL database for carrying out the tests.
This database is cleaned up after testing so it can be reused.

The test database needs to be created first with the CREATE DATABASE command;
this is not done in the testing environment.
"""

import unittest
import os
import biom
import psycopg2
import networkx as nx
from masq.scripts.sq4biom import BiomConnection
from masq.scripts.io import IoConnection
from masq.scripts.netstats import SetConnection, _convert_network, extract_sets


__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

testraw = """{
     "id":  "test",
     "format": "Biological Observation Matrix 1.0.0-dev",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteoba\
cteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriac\
eae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobact\
eria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichosp\
ermum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchae\
ota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "\
g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaer\
obium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
        ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 6],
     "data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
            ]
    }
"""

testbiom = biom.parse.parse_biom_table(testraw)
testdict = dict.fromkeys(testbiom._observation_ids)

# make toy network
g1 = nx.Graph()
nodes = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]
g1.add_nodes_from(nodes)
g1.add_edges_from([("GG_OTU_1", "GG_OTU_2"),
                  ("GG_OTU_2", "GG_OTU_5"), ("GG_OTU_3", "GG_OTU_4")])
g1["GG_OTU_1"]["GG_OTU_2"]['weight'] = 1.0
g1["GG_OTU_2"]["GG_OTU_5"]['weight'] = 1.0
g1["GG_OTU_3"]["GG_OTU_4"]['weight'] = -1.0

g2 = g1.copy(as_view = False)
g2.remove_edge("GG_OTU_3", "GG_OTU_4")
g2.add_edge("GG_OTU_3", "GG_OTU_5", weight=-1.0)

# check for edge weight
g2.edges[("GG_OTU_1", "GG_OTU_2")]['weight'] = -1.0
g1.add_edge("GG_OTU_5", "GG_OTU_3", weight=1.0)


class TestNetstats(unittest.TestCase):
    """
    Tests netstats methods.
    Warning: most of these functions are to interact with a local database named test.
    Therefore, the presence of the necessary local files is a prerequisite.
    """
    @classmethod
    def setUpClass(cls):
        # The class setup creates a config file
        # this config file refers to the local test database
        # if your test database has different config, change here
        config = "[postgresql]\n" \
                 "host=localhost\n" \
                 "database=test\n" \
                 "user=test\n" \
                 "password=test\n"
        file = open("database.ini", "w")
        file.write(config)
        file.close()
        # clear database before usage
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        tables = ['bioms', 'sample', 'networks', 'taxonomy', 'counts', 'edges', 'meta']
        for tab in tables:
            try:
                cur.execute(("DROP TABLE " + tab + ";"))
            except psycopg2.Error:
                pass
        conn.commit()
        cur.close()
        conn.close()

    @classmethod
    def tearDownClass(cls):
        os.remove("database.ini")
        # clear database after usage
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        tables = ['bioms', 'sample', 'networks', 'taxonomy', 'counts', 'edges', 'meta']
        for tab in tables:
            try:
                cur.execute(("DROP TABLE " + tab + ";"))
            except psycopg2.Error:
                pass
        conn.commit()
        cur.close()
        conn.close()

    def test_extract_sets(self):
        """
        Tests if the import_networks function reads the correct database file,
        and imports the network file, by querying the networks table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        extract_sets(set='intersection', path=os.getcwd(), weight=False)
        file = nx.read_graphml("intersection.graphml")
        conn_object.delete_tables()
        os.remove("intersection.graphml")
        self.assertEqual(len(file.edges), 3)

    def test_get_intersection(self):
        """
        Tests if the import_network function reads the correct database file,
        and imports the network file, by querying the edges table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        intersection_set = conn_object.get_intersection(networks=["g1", "g2"], number=2)
        conn_object.delete_tables()
        self.assertCountEqual(("GG_OTU_2", "GG_OTU_5"), list(intersection_set.edges)[0])

    def test_get_intersection_3(self):
        """
        Tests if the intersection returns an empty network
        when the intersection is too large.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        intersection_set = conn_object.get_intersection(networks=["g1", "g2"], number=3)
        conn_object.delete_tables()
        self.assertEqual(len(intersection_set.edges), 0)

    def test_get_intersection_weight(self):
        """
        Tests if the intersection returns the correct number of edges
        when the weight parameter is changed.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        intersection_set = conn_object.get_intersection(networks=["g1", "g2"], number=2, weight=False)
        conn_object.delete_tables()
        self.assertEqual(len(intersection_set.edges), 5)

    def test_get_difference(self):
        """
        Tests if the difference is returned with edges that have different weights.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        difference_set = conn_object.get_difference(networks=["g1", "g2"])
        conn_object.delete_tables()
        self.assertEqual(len(difference_set.edges), 3)

    def test_get_difference_weight(self):
        """
        Tests if the difference is returned with only edges that have the same weights.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        difference_set = conn_object.get_difference(networks=["g1", "g2"], weight=False)
        conn_object.delete_tables()
        self.assertEqual(len(difference_set.edges), 1)

    def test_get_union(self):
        """
        Tests if the union of the networks is returned.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        union_set = conn_object.get_union(networks=["g1", "g2"])
        conn_object.delete_tables()
        self.assertTrue(len(union_set.edges), 4)

    def test_aggr_networks(self):
        """
        Tests if network names are aggregated to an edge property
        named 'source'.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        difference_set = conn_object.get_difference(networks=["g1", "g2"])
        conn_object.delete_tables()
        self.assertEqual(difference_set.edges[('GG_OTU_1', 'GG_OTU_2')]['source'], 'g1')

    def test_aggr_weight(self):
        """
        Tests if edge weights are aggregated to an edge property
        named 'weight'.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        conn_object = SetConnection()
        intersection_set = conn_object.get_intersection(networks=["g1", "g2"], number=2, weight=False)
        conn_object.delete_tables()
        self.assertEqual(intersection_set.edges[('GG_OTU_1', 'GG_OTU_2')]['weight'], '1,-1')

    def test_SetConnection(self):
        """
        Tests if the network file is added
        to the database.
        :return:
        """
        test_config = {'host': 'localhost',
                       'database': 'test',
                       'user': 'test',
                       'password': 'test'}
        conn = SetConnection()
        read_config = conn.config
        self.assertTrue(test_config == read_config)

    def test_convert_network(self):
        """
        Tests whether the set result is converted to a Networkx graph with 1 edge.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g1, name='g1', study='test')
        conn_object.add_network(network=g2, name='g2', study='test')
        set_query = "SELECT string_agg(networkID::varchar, ',') AS networks, " \
                    "source, target, SIGN(weight), COUNT(*) " \
                    "FROM edges " \
                    "WHERE networkID IN " + str(('g1', 'g2')) + \
                    " GROUP BY source, target, SIGN(weight) " \
                    "HAVING COUNT(*) > " + str(1)
        set_result = self.query(set_query, fetch=True)
        graph = _convert_network(set_result)
        conn_object.delete_tables()
        self.assertEqual(len(graph.edges), 1)


if __name__ == '__main__':
    unittest.main()
