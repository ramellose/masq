"""
This file contains functions for testing functions
and the MetaConnection class in the metastats.py file.

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
from masq.scripts.io import IoConnection
from masq.scripts.sq4biom import BiomConnection
from masq.scripts.metastats import start_metastats, MetaConnection
from pandas.core.common import flatten


__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

testraw = """{
     "id":null,
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
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Punk", "\
g_Anthrax", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Punk", "g__NOFX", "s__"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}, 
        {"id":"GG_OTU_6", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Punk", "g__Misfits", "s__"]}}\
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "pH":"2.0",
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "pH":"1.8",       
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "pH":"2.3",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample5", "metadata":{
                                "pH":"2.0",        
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample6", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample7", "metadata":{
                                "pH":"1.9",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample8", "metadata":{
                                "pH":"1.9",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},    
        {"id":"Sample9", "metadata":{
                                "pH":"1.8",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample10", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},                                    
        {"id":"Sample11", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample12", "metadata":{
                                "pH":"6.9",        
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample13", "metadata":{
                                "pH":"7.1",        
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample14", "metadata":{
                                "pH":"7.0",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample15", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample16", "metadata":{
                                "pH":"6.9",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample17", "metadata":{
                                "pH":"6.7",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},   
        {"id":"Sample18", "metadata":{
                                "pH":"7.2",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample19", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},                                                                                                                         
        {"id":"Sample20", "metadata":{
                                "pH":"7.0",        
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
        ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 20],
     "data":[[0,10,5],
             [0,11,5],
             [0,12,6],
             [0,13,5],
             [0,14,5],
             [0,15,5],
             [0,16,6],
             [0,17,5],
             [0,18,5],
             [0,19,6],
             [0,9,6],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,8,5],
             [1,10,1],
             [1,11,2],
             [1,2,3],
             [1,14,5],
             [1,17,1],
             [1,12,2],
             [1,19,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [2,6,1],
             [2,8,4],
             [2,10,2],
             [2,14,4],
             [2,16,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [3,7,2],
             [3,12,1],
             [3,15,2],
             [3,7,1],
             [3,10,1],
             [3,11,1],
             [4,1,1],
             [4,2,1],
             [4,4,1],
             [4,14,1],
             [4,6,1]
            ]
    }
"""

testbiom = biom.parse.parse_biom_table(testraw)
testdict = dict.fromkeys(testbiom._observation_ids)

# make toy network
g = nx.Graph()
nodes = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]
g.add_nodes_from(nodes)
g.add_edges_from([("GG_OTU_1", "GG_OTU_3"),
                  ("GG_OTU_2", "GG_OTU_5"), ("GG_OTU_4", "GG_OTU_5"),
                  ("GG_OTU_2", "GG_OTU_6")])
g["GG_OTU_1"]["GG_OTU_3"]['weight'] = 1.0
g["GG_OTU_2"]["GG_OTU_5"]['weight'] = 1.0
g["GG_OTU_4"]["GG_OTU_5"]['weight'] = -1.0
g["GG_OTU_2"]["GG_OTU_6"]['weight'] = -1.0


class TestMeta(unittest.TestCase):
    """
    Tests IO methods.
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

    def test_start_metastats(self):
        """
        Tests if the stat_networks function reads the correct database file,
        and agglomerates the network file, by querying the networks table.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        start_metastats(level='Family', weight=False, networks=['g'])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT node_num "
                    "FROM networks WHERE networks.networkid = %s"
                    "LIMIT 1;", vars=('Genus_g',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 5)

    def test_agglomerate_networks(self):
        """
        Tests if the agglomerate_networks function agglomerates the test file
        so that only 3 nodes and 3 edges remain.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        conn_object.agglomerate_networks(level='Family', weight=False, networks=['g'])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT node_num "
                    "FROM networks WHERE networks.networkid = %s"
                    "LIMIT 1;", vars=('Family_g',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 3)

    def test_get_pairlist(self):
        """
        Tests if the pair list is returned correctly,
        with a pair being 2 edges with matching partners at the specified taxonomic level.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        pair = conn_object.get_pairlist(level='Family', weight=False, network='g')
        conn_object.delete_tables()
        self.assertCountEqual(pair[0], ['GG_OTU_1', 'GG_OTU_4'])

    def test_get_pairlist_weight_noresults(self):
        """
        Tests if the pair is not returned when weight is set to True.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        pair = conn_object.get_pairlist(level='Family', weight=True, network='g')
        conn_object.delete_tables()
        self.assertEqual(len(pair), 0)

    def test_get_pairlist_weight_esults(self):
        """
        Tests if the pair is returned when weight is set to True and edge weights set to the same. .
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        f = g.copy(as_view=False)
        f.edges[("GG_OTU_4","GG_OTU_5")]['weight'] = 1.0
        conn_object.add_network(network=f, name='g', study='test')
        conn_object = MetaConnection()
        pair = conn_object.get_pairlist(level='Family', weight=True, network='g')
        conn_object.delete_tables()
        self.assertEqual(len(pair), 3)

    def test_get_taxlist(self):
        """
        Tests if the matching taxa with same taxonomy at Family level are returned.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        pair = conn_object.get_taxlist(level='Family', network='g')
        conn_object.delete_tables()
        self.assertCountEqual(pair, ['GG_OTU_3', 'GG_OTU_4', 'GG_OTU_6'])

    def test_MetaConnection(self):
        """
        Tests if the IoConnection can be successfully initialized,
        with the correct config parameters.
        :return:
        """
        test_config = {'host': 'localhost',
                       'database': 'test',
                       'user': 'test',
                       'password': 'test'}
        conn = MetaConnection()
        read_config = conn.config
        self.assertTrue(test_config == read_config)

    def test_copy_network(self):
        """
        Tests if the network file is added
        to the database.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        conn_object.copy_network(source_network='g', new_network='f')
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT networkID "
                    "FROM networks WHERE networks.networkid = %s"
                    "LIMIT 1;", vars=('f',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'f')

    def test_agglomerate_taxa(self):
        """
        Tests whether two taxa are merged into one node.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        conn_object.agglomerate_taxa(['GG_OTU_1', 'GG_OTU_2'], network='g', level='Class')
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT source, target "
                    "FROM edges WHERE edges.networkid = %s;", vars=('g',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        result = list(flatten(result))
        removed = {"GG_OTU_1", "GG_OTU_2"}
        self.assertEqual(len(removed.intersection(result)), 0)

    def test_agglomerate_pair(self):
        """
        Tests whether two edges are merged into a single edge,
        so that OTU 1 and 3 are removed.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        conn_object.agglomerate_pair((["GG_OTU_1", "GG_OTU_2"],
                                      ["GG_OTU_3", "GG_OTU_5"], 1),
                                     network='g', level='Class')
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT source, target "
                    "FROM edges WHERE edges.networkid = %s;", vars=('g',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        removed = {"GG_OTU_1", "GG_OTU_3"}
        self.assertEqual(len(removed.intersection(result)), 0)

    def test_create_agglom(self):
        """
        Tests whether a new agglomerated node is created with intact taxonomy.
        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        uid = conn_object.create_agglom(parent="GG_OTU_1", level="Family")
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT Family "
                    "FROM taxonomy WHERE taxonomy.taxon = %s;", vars=(uid,))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 'f__Enterobacteriaceae')

    def test_agglomerate_network_phylumloop(self):
        """
        Tests whether the network agglomeration also works up to Phylum level.
        Note: need to check selfloops.

        :return:
        """
        conn_object = BiomConnection()
        conn_object.create_tables()
        conn_object.add_biom(testbiom, 'test')
        conn_object = IoConnection()
        conn_object.add_network(network=g, name='g', study='test')
        conn_object = MetaConnection()
        conn_object.agglomerate_networks(level='Phylum', weight=False, networks=['g'])
        conn = psycopg2.connect(**{"host": "localhost",
                                   "database": "test",
                                   "user": "test",
                                   "password": "test"})
        cur = conn.cursor()
        cur.execute("SELECT node_num "
                    "FROM networks WHERE networks.networkid = %s"
                    "LIMIT 1;", vars=('Phylum_g',))
        result = cur.fetchall()
        cur.close()
        conn.close()
        conn_object.delete_tables()
        self.assertEqual(result[0][0], 3)


if __name__ == '__main__':
    unittest.main()
