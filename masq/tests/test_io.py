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
from biom.cli.util import write_biom_table
from masq.scripts.io import import_networks, IoConnection, _read_network_extension
from masq.scripts.sq4biom import BiomConnection


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
g = nx.Graph()
nodes = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]
g.add_nodes_from(nodes)
g.add_edges_from([("GG_OTU_1", "GG_OTU_2"),
                  ("GG_OTU_2", "GG_OTU_5"), ("GG_OTU_3", "GG_OTU_4")])
g["GG_OTU_1"]["GG_OTU_2"]['weight'] = 1.0
g["GG_OTU_2"]["GG_OTU_5"]['weight'] = 1.0
g["GG_OTU_3"]["GG_OTU_4"]['weight'] = -1.0


class TestIo(unittest.TestCase):
    """
    Tests utils methods.
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

    def test_import_networks(self):
        """
        Tests if the import_networks function reads the correct database file,
        and imports the network file, by querying the networks table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(('test', 2, 5))
        nx.write_graphml(g, path="test.graphml")
        import_networks("test.graphml")
        os.remove("test.graphml")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM networks "
                    "LIMIT 1;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'test')

    def test_import_network_edge(self):
        """
        Tests if the import_network function reads the correct database file,
        and imports the network file, by querying the edges table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(('test', 2, 5))
        nx.write_graphml(g, path="test.graphml")
        import_networks("test.graphml")
        os.remove("test.graphml")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT source, target "
                    "FROM edges "
                    "LIMIT 2;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertCountEqual(result,
                              [("GG_OTU_1", "GG_OTU_2"),
                                ("GG_OTU_2", "GG_OTU_5")])

    def test_import_network_mapping(self):
        """
        Tests if the IoConnection uses the mapping dict.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(('banana', 2, 5))
        nx.write_graphml(g, path="test.graphml")
        import_networks("test.graphml", mapping={"test": "banana"})
        os.remove("test.graphml")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM networks "
                    "LIMIT 1;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'banana')

    def test_import_network_source(self):
        """
        Tests if the IoConnection uses the source dict.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(('banana', 2, 5))
        nx.write_graphml(g, path="test.graphml")
        import_networks("test.graphml", sources={"test": "banana"})
        os.remove("test.graphml")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM networks "
                    "LIMIT 1;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'test')

    def test_IoConnection(self):
        """
        Tests if the IoConnection can be successfully initialized,
        with the correct config parameters.
        :return:
        """
        test_config = {'host': 'localhost',
                       'database': 'test',
                       'user': 'test',
                       'password': 'test'}
        conn = IoConnection()
        read_config = conn.config
        self.assertTrue(test_config == read_config)

    def test_add_network(self):
        """
        Tests if the network file is added
        to the database.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'banana')
        conn_object = IoConnection()
        conn_object.add_network(g, 'banana', 'banana')
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from edges;")
        result = cur.fetchall()
        # 3 edges in g
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(len(result), 3)

    def test_add_network_node(self):
        """
        Tests whether a row is added to the networks table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'potato')
        conn_object = IoConnection()
        conn_object.add_network_node(values=('potato', 'potato', 5, 20))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from networks;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0], ('potato', 'potato', 5, 20))

    def test_add_edge(self):
        """
        Tests whether a row is added to the edge table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'banana')
        conn_object = IoConnection()
        conn_object.add_edge(values=("banana", "GG_OTU_1", "GG_OTU_2", 0.3))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from edges;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0], ("banana", "GG_OTU_1", "GG_OTU_2", 0.3))

    def test_add_edges(self):
        """
        Tests whether new rows are added to the edges table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'banana')
        conn_object = IoConnection()
        conn_object.add_edge(values=[("banana", "GG_OTU_1", "GG_OTU_2", 0.3),
                                     ("banana", "GG_OTU_3", "GG_OTU_5", -0.8)])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from edges;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(len(result), 2)

    def test_edge_error(self):
        """
        Tests if the query correctly reports an error
        when network data is not present
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'banana')
        conn_object = IoConnection()
        self.assertRaises(psycopg2.Error, conn_object.add_edge(
            values=('apple', 'Streptococcus', 'Sample1', 0.5)),
        )
        # long format table of 5 * 6 observations
        conn_object.delete_tables()

    def test_read_network_extension(self):
        """
        Tests whether the network can be read from a file.
        :return:
        """
        nx.write_graphml(g, path="test.graphml")
        graph = _read_network_extension("test.graphml")
        os.remove("test.graphml")
        self.assertEqual(len(graph.edges), 3)


if __name__ == '__main__':
    unittest.main()
