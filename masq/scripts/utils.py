"""
This file contains a parent class for setting up a connection
to a sqlite3 database,
and populating this database with tables that can contain count data,
network data and (taxon and sample) metadata. """

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import psycopg2
from configparser import ConfigParser
import sys
import os
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


def setup_database(config='database.ini', create=True,
                   host=None, database=None,
                   username=None, password=None):
    """
    Connects to PostgreSQL database,
    either via provided args (priority)
    or the database.ini config file.
    Depending on the arguments,
    tables for adding BIOM files and networks
    are either added or deleted.

    :param config: Location of file with database parameters.
    :param create: If true, tables are created. Otherwise, deleted.
    :param host: Database address.
    :param database: Name of PostgreSQL database.
    :param username: Username for PostgreSQL database.
    :param password: Password of PostgreSQL database.
    :return:
    """
    conn = ParentConnection(config, host, database, username, password)
    if create:
        conn.create_tables()
    else:
        conn.delete_tables()


class ParentConnection:
    def __init__(self, config='database.ini',
                 host=None, database=None,
                 username=None, password=None):
        """
        Initialzes a driver for accessing the PostgreSQL database.
        The driver can either read in a config file,
        or accept parameters directly.

        Config adapted from: https://www.postgresqltutorial.com/postgresql-python/connect/

        :param config: Location of file with database parameters.
        :param host: Database address.
        :param database: Name of PostgreSQL database.
        :param username: Username for PostgreSQL database.
        :param password: Password of PostgreSQL database.
        """
        conn = None
        self.config = {'host': None,
                       'database': None,
                       'user': None,
                       'password': None}
        new_params = {'host': host,
                       'database': database,
                       'user': username,
                       'password': password}
        if config:
            parser = ConfigParser()
            try:
                if os.path.isfile(config):
                    parser.read(config)
                else:
                    parser.read(os.getcwd() + "\\" + config)
            except FileNotFoundError:
                logger.warning("Could not read database config file.")
                exit()
            params = parser.items('postgresql')
            for param in params:
                self.config[param[0]] = param[1]
        for d in new_params:
            if new_params[d]:
                self.config[d] = new_params[d]
        try:
            logger.info("Connecting to the PostgreSQL database...")
            conn = psycopg2.connect(**self.config)
            cur = conn.cursor()
            cur.execute("SELECT version()")
            db_version = cur.fetchone()
            logger.info("Running PostgreSQL version " + db_version[0])
            cur.close()
        except psycopg2.DatabaseError as e:
            logger.warning(e)
        conn.close()

    def query(self, query, fetch=False):
        """
        Accepts a sqlite query and provides the results (??)
        The connection is opened and closed within this query,
        so the connection is always closed safely.
        :param query: String containing sqlite query
        :param fetch: If set to true, fetches output
        :return:
        """
        conn = psycopg2.connect(**self.config)
        results = None
        try:
            c = conn.cursor()
            c.execute(query)
            if fetch:
                results = c.fetchall()
        except psycopg2.Error as e:
            logger.error(e)
        conn.commit()
        c.close()
        conn.close()
        return results

    def value_query(self, query, values, fetch=False):
        """
        If given a standard query with question marks,
        the values are used to replace the question marks in the query.
        This query is an easy way to safely insert rows.
        If the values are a tuple, a single query is executed.
        If values are a list of tuples, many queries are executed.

        :param query: String containing sqlite query with question marks
        :param values: Iterable of values
        :param fetch: If set to true, fetches output
        :return: Last row ID
        """
        conn = psycopg2.connect(**self.config)
        results = None
        if type(values) == tuple:
            try:
                c = conn.cursor()
                c.execute(query, values)
            except psycopg2.Error as e:
                logger.error(e)
        elif type(values) == list:
            try:
                c = conn.cursor()
                c.executemany(query, values)
            except psycopg2.Error as e:
                logger.error(e)
        else:
            logger.warning("Values are not a tuple or list, so no query was executed.")
        if fetch:
            results = c.fetchall()
        conn.commit()
        c.close()
        conn.close()
        return results

    def create_tables(self):
        """
        Adds the default data entities.
        :return:
        """
        biom_query = "CREATE TABLE IF NOT EXISTS bioms (" \
                     "studyID varchar PRIMARY KEY," \
                     "tax_num int," \
                     "sample_num int" \
                     ");"
        network_query = "CREATE TABLE IF NOT EXISTS networks (" \
                        "networkID varchar PRIMARY KEY," \
                        "studyID varchar," \
                        "node_num int," \
                        "edge_num int," \
                        "FOREIGN KEY (studyID) REFERENCES bioms(studyID) ON DELETE CASCADE" \
                        ");"
        tax_query = "CREATE TABLE IF NOT EXISTS taxonomy (" \
                    "taxon varchar PRIMARY KEY," \
                    "studyID varchar," \
                    "Kingdom varchar," \
                    "Phylum varchar," \
                    "Class varchar," \
                    '"Order" varchar,' \
                    "Family varchar," \
                    "Genus varchar," \
                    "Species varchar," \
                    "FOREIGN KEY (studyID) REFERENCES bioms(studyID) ON DELETE CASCADE" \
                    ");"
        edge_query = "CREATE TABLE IF NOT EXISTS edges (" \
                     "networkID varchar," \
                     "source varchar NOT NULL," \
                     "target varchar NOT NULL," \
                     "weight float," \
                     "FOREIGN KEY (networkID) REFERENCES networks(networkID) ON DELETE CASCADE" \
                     ");"
        # in the metadata table,
        # we actually need 2 tables:
        # one with unique meta id + study id,
        # other with property values.
        # leave one column empty if value is not appliccable
        # the metadata table needs to be 'long' format;
        # this to tackle issues with different tables
        # having different columns.
        meta_id_query = "CREATE TABLE IF NOT EXISTS samples (" \
                        "sampleID varchar PRIMARY KEY," \
                        "studyID varchar," \
                        "FOREIGN KEY (studyID) REFERENCES bioms(studyID) ON DELETE CASCADE" \
                        ");"
        meta_query = "CREATE TABLE IF NOT EXISTS meta (" \
                     "sampleID varchar," \
                     "studyID varchar," \
                     "property varchar," \
                     "textvalue varchar," \
                     "numvalue float," \
                     "FOREIGN KEY (studyID) REFERENCES bioms(studyID) ON DELETE CASCADE," \
                     "FOREIGN KEY (sampleID) REFERENCES samples(sampleID)" \
                     ");"
        counts_query = "CREATE TABLE IF NOT EXISTS counts (" \
                       "studyID varchar," \
                       "taxon varchar," \
                       "sampleID varchar," \
                       "count float," \
                       "FOREIGN KEY (SampleID) REFERENCES samples(sampleID)," \
                       "FOREIGN KEY (Taxon) REFERENCES taxonomy(Taxon)," \
                       "FOREIGN KEY (studyID) REFERENCES bioms(studyID) ON DELETE CASCADE" \
                       ");"
        for query in [biom_query, network_query, tax_query, edge_query,
                      meta_id_query, meta_query, counts_query]:
            try:
                self.query(query)
            except psycopg2.Error as e:
                logger.warning(e)

    def delete_tables(self):
        """
        Deletes the tables created by the create_tables function.
        :return:
        """
        self.query("DROP TABLE counts;")
        self.query("DROP TABLE edges;")
        self.query("DROP TABLE meta;")
        self.query("DROP TABLE samples;")
        self.query("DROP TABLE taxonomy;")
        self.query("DROP TABLE networks CASCADE;")
        self.query("DROP TABLE bioms CASCADE;")

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