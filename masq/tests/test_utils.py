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
from masq.scripts.utils import ParentConnection


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


class TestUtils(unittest.TestCase):
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

    def test_ParentConnection(self):
        """
        Tests if the ParentConnection can be successfully initialized,
        with the correct config parameters.
        :return:
        """
        test_config = {'host': 'localhost',
                       'database': 'test',
                       'user': 'test',
                       'password': 'test'}
        conn = ParentConnection()
        read_config = conn.config
        self.assertTrue(test_config == read_config)

    def test_create_tables(self):
        """
        Tests if the requested tables are created.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from information_schema.tables "
                    "WHERE table_schema = 'public' "
                    "AND table_type = 'BASE TABLE';")
        result = cur.fetchall()
        result = [x[2] for x in result]
        self.assertCountEqual(result,
                             ['bioms', 'counts', 'networks',
                              'taxonomy', 'edges', 'samples', 'meta']
                             )
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_delete_tables(self):
        """
        Tests whether there are no tables in the database
        if the object methods add and then delete these.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        conn_object.delete_tables()
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from information_schema.tables "
                    "WHERE table_schema = 'public' "
                    "AND table_type = 'BASE TABLE';")
        result = cur.fetchall()
        self.assertEqual(len(result), 0)
        cur.close()
        conn.close()

    def test_query(self):
        """
        Tests whether the query accurately deletes the meta table.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        conn_object.query("DROP TABLE meta")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from information_schema.tables "
                    "WHERE table_schema = 'public' "
                    "AND table_type = 'BASE TABLE';")
        result = cur.fetchall()
        result = [x[2] for x in result]
        self.assertFalse('meta' in result)
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_value_query(self):
        """
        Tests whether the query adds a row to the bioms table,
        by first adding the row and then getting the studyID.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        biom_query = "INSERT INTO bioms (studyID,tax_num,sample_num) " \
                       "VALUES (%s,%s,%s)"
        conn_object.value_query(biom_query, values=("test", 300, 200))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM bioms "
                    "LIMIT 1;")
        result = cur.fetchall()
        self.assertEqual(result[0][0], 'test')
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_values_query(self):
        """
        Tests whether the query adds several rows to the bioms table,
        by first adding the rows and then getting the studyID.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        biom_query = "INSERT INTO bioms (studyID,tax_num,sample_num) " \
                       "VALUES (%s,%s,%s)"
        conn_object.value_query(biom_query, values=[("test", 300, 200),
                                                    ("test2", 400, 1500)])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID,tax_num,sample_num "
                    "FROM bioms;")
        result = cur.fetchall()
        self.assertListEqual([('test', 300, 200),
                              ('test2', 400, 1500)],
                             result)
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_value_query_error(self):
        """
        Tests if the value query correctly reports an error
        when a list is given instead of a tuple,
        so item lengths no longer match and the SQL query crashes.
        :return:
        """
        conn_object = ParentConnection()
        conn_object.create_tables()
        biom_query = "INSERT INTO bioms (studyID,tax_num,sample_num) " \
                     "VALUES (%s,%s,%s)"
        self.assertRaises(TypeError, conn_object.value_query,
                          biom_query,
                          ["test", 300, 200])
        conn_object.delete_tables()


if __name__ == '__main__':
    unittest.main()
