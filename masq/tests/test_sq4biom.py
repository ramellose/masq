"""
This file contains functions for testing functions
and the BiomConnection class in the sq4biom.py file.

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
from masq.scripts.sq4biom import BiomConnection, import_biom


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


class TestBiom(unittest.TestCase):
    """
    Tests sq4biom methods.
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

    def test_import_biom(self):
        """
        Tests if the import_biom function reads the correct database file,
        and imports the biom file, by querying the bioms table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        write_biom_table(testbiom, fmt='hdf5', filepath="test.biom")
        import_biom("test.biom", mapping=None)
        os.remove("test.biom")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM bioms "
                    "LIMIT 1;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'test')

    def test_import_biom_obs(self):
        """
        Tests if the import_biom function reads the correct database file,
        and imports the biom file, by queriying the counts table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        write_biom_table(testbiom, fmt='hdf5', filepath="test.biom")
        import_biom("test.biom", mapping=None)
        os.remove("test.biom")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT sampleid "
                    "FROM counts "
                    "LIMIT 5;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        result = [x[0] for x in result]
        self.assertCountEqual(result,
                              ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'])

    def test_import_biom_mapping(self):
        """
        Tests if the BiomConnection uses the mapping dict.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        write_biom_table(testbiom, fmt='hdf5', filepath="test.biom")
        import_biom("test.biom", mapping={'test': 'banana'})
        os.remove("test.biom")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT studyID "
                    "FROM bioms "
                    "LIMIT 1;")
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'banana')

    def test_BiomConnection(self):
        """
        Tests if the BiomConnection can be successfully initialized,
        with the correct config parameters.
        :return:
        """
        test_config = {'host': 'localhost',
                       'database': 'test',
                       'user': 'test',
                       'password': 'test'}
        conn = BiomConnection()
        read_config = conn.config
        self.assertTrue(test_config == read_config)

    def test_add_biom(self):
        """
        Tests if the BIOM file is imported by
        querying the count data.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'banana')
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from counts;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(len(result), 30)
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_summary(self):
        """
        Tests whether a row is added to the bioms table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('potato', 5, 20))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from bioms;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0], ('potato', 5, 20))
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_taxon(self):
        """
        Tests whether a row is added to the taxonomy table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_taxon(values=('Listeria', 'banana', 'Bacteria',
                                      'Firmicutes', 'Bacilli',
                                      'Bacillales', 'Listeriaceae',
                                      'Listeria', 'monocytogenes'))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from taxonomy;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0], ('Listeria', 'banana',
                                     'Bacteria',
                                     'Firmicutes', 'Bacilli',
                                     'Bacillales', 'Listeriaceae',
                                     'Listeria', 'monocytogenes'))
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_taxon_missingvalues(self):
        """
        Tests whether a row is added to the taxonomy table,
        even when some assignments are unknown.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_taxon(values=('Listeria', 'banana', 'Bacteria',
                                      'Firmicutes', 'Bacilli',
                                      None, None, None, None))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT * from taxonomy;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0], ('Listeria', 'banana',
                                     'Bacteria',
                                     'Firmicutes', 'Bacilli',
                                     None, None, None, None))
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_sample(self):
        """
        Tests whether a row is added to the sample table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_sample(values=('Sample1', 'banana'))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT sampleid from samples;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0][0], 'Sample1')
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_meta(self):
        """
        Tests whether a row is added to the meta table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_sample(values=('Sample1', 'banana'))
        conn_object.add_meta(values=('Sample1', 'banana',
                                     'colour', 'yellow', None))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT property from meta;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0][0], 'colour')
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_observation(self):
        """
        Tests whether a row is added to the counts table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_sample(values=('Sample1', 'banana'))
        conn_object.add_taxon(values=('Listeria', 'banana', 'Bacteria',
                                      'Firmicutes', 'Bacilli',
                                      'Bacillales', 'Listeriaceae',
                                      'Listeria', 'monocytogenes'))
        conn_object.add_observation(values=('banana', 'Listeria',
                                            'Sample1', 0.5))
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT count from counts;")
        result = cur.fetchall()
        # long format table of 5 * 6 observations
        self.assertEqual(result[0][0], 0.5)
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_add_observations(self):
        """
        Tests whether multiple rows are added to the counts table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_sample(values=('Sample1', 'banana'))
        conn_object.add_taxon(values=[('Listeria', 'banana', 'Bacteria',
                                      'Firmicutes', 'Bacilli',
                                      'Bacillales', 'Listeriaceae',
                                      'Listeria', 'monocytogenes'),
                                      ('Listeria2', 'banana', 'Bacteria',
                                       'Firmicutes', 'Bacilli',
                                       'Bacillales', 'Listeriaceae',
                                       'Listeria', 'monocytogenes')
                                      ])
        conn_object.add_observation(values=[('banana', 'Listeria',
                                            'Sample1', 0.5),
                                            ('banana', 'Listeria2',
                                             'Sample1', 0.75)
                                            ])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT count from counts;")
        result = cur.fetchall()
        self.assertEqual(len(result), 2)
        cur.close()
        conn.close()
        conn_object.delete_tables()

    def test_observation_error(self):
        """
        Tests if the query correctly reports an error
        when a taxon is not present in the taxonomy table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_summary(values=('banana', 1, 2))
        conn_object.add_sample(values=('Sample1', 'banana'))
        conn_object.add_taxon(values=('Listeria', 'banana', 'Bacteria',
                                      'Firmicutes', 'Bacilli',
                                      'Bacillales', 'Listeriaceae',
                                      'Listeria', 'monocytogenes'))
        self.assertRaises(psycopg2.Error, conn_object.add_observation(
            values=('banana', 'Streptococcus', 'Sample1', 0.5)),
        )
        # long format table of 5 * 6 observations
        conn_object.delete_tables()


if __name__ == '__main__':
    unittest.main()
