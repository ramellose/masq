"""
This file contains a function for reading in command-line arguments
so BIOM files or tab-delimited files can be read.
The BIOM format is a standardized biological format
that is commonly used to contain biological data.
Tab-delimited files should be supplied with the BIOM-format specified headers (# prefix).

The software can operate in two manners:
import all BIOM files in a folder,
or import separate BIOM files / tab-delimited files

The file also defines a class for a Neo4j driver.
Given a running database, this driver can upload and delete experiments in the database.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import os
import sys
import biom
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


def import_biom(location, mapping=None,
                config='database.ini',
                host=None, database=None,
                username=None, password=None):
    """
    Can read a single BIOM file or all BIOM files in a folder.
    These are then imported into the sqlite database.
    If no mapping is supplied,
    the filename is used to match the BIOM file to the network.

    :param location: Location of BIOM file(s)
    :param mapping: Dictionary with BIOM filenames as keys, new names as values
    :param config: Location of file with database parameters.
    :param host: Database address.
    :param database: Name of PostgreSQL database.
    :param username: Username for PostgreSQL database.
    :param password: Password of PostgreSQL database.
    :return:
    """
    conn = BiomConnection(config, host, database, username, password)
    if os.path.isdir(location):
        for y in os.listdir(location):
            biomtab = biom.load_tabe(location + '/' + y)
            name = y.split(".")[0]
            if mapping:
                name = mapping[name]
            conn.add_biom(biomtab, name)
    else:
        biomtab = biom.load_tabe(location)
        name = location.split('/')[-1]
        name = name.split('\\')[-1]
        name = name.split(".")[0]
        if mapping:
            name = mapping[name]
        conn.add_network(biomtab, name)


class BiomConnection(ParentConnection):
    """
    Initializes a connection to the sqlite3 database.
    This connection object contains methods for converting NetworkX
    objects to rows in the summary network table
    and rows in the edges table.
    """
    # inherits init from parent
    def add_biom(self, biomfile, name):
        """
        Takes a networkx object and writes this to the sqlite3 database.
        :param biomfile: BIOM table
        :param name: BIOM table name
        :return:
        """
        values = (name, biomfile.shape[0], biomfile.shape[1])
        self.add_summary(values)
        taxa = biomfile.ids(axis='observation')
        samples = biomfile.ids(axis='sample')
        taxonomy_values = list()
        sample_values = list()
        meta_values = list()
        obs_values = list()
        for sample in samples:
            values = list()
            sample_values.append((sample, name))
        for tax in taxa:
            taxonomy = biomfile.metadata(id=tax, axis='observation')['taxonomy']
            values = [tax, name]
            values.extend(taxonomy)
            # last values may be removed if taxonomy is unavailable
            while len(values) < 9:
                values.append(None)
            taxonomy_values.append(tuple(values))
            data = biomfile.data(id=tax, axis='observation')
            for sample in samples:
                sample_data = biomfile.metadata(id=sample, axis='sample')
                values = list()
                if sample_data:
                    for property in sample_data:
                        value = [sample]
                        value.append(name)
                        value.append(property)
                        if type(sample_data[property]) == float or type(sample_data[property]) == int:
                            value.append(None)
                            value.append(sample_data[property])
                        else:
                            value.append(sample_data[property])
                            value.append(None)
                        values.append(tuple(value))
                    meta_values.extend(values)
                count_index = biomfile.index(sample, axis='sample')
                count = data[count_index]
                values = name, tax, sample, count
                obs_values.append(values)
        self.add_taxon(taxonomy_values)
        self.add_sample(sample_values)
        self.add_meta(meta_values)
        self.add_observation(obs_values)
        logger.info("Uploaded BIOM data for " + name +".\n")

    def add_summary(self, values):
        """
        Adds rows to the bioms table in the PostgreSQL database.

        :param values: List of tuples for bioms table
        :return:
        """
        biom_query = "INSERT INTO bioms (studyID,tax_num,sample_num) " \
                       "VALUES (%s,%s,%s)"
        self.value_query(biom_query, values)

    def add_taxon(self, values):
        """
        Adds taxonomy rows to the taxonomy table in the PostgreSQL database.

        :param values: List of tuples for taxonomy table
        :return:
        """
        tax_query = 'INSERT INTO taxonomy (taxon,studyID,Kingdom,Phylum,Class,"Order",Family,Genus,Species) ' \
                    'VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)'
        self.value_query(tax_query, values)

    def add_sample(self, values):
        """
        Adds rows to the sample table in the PostgreSQL database.

        :param values: List of tuples for sample table
        :return:
        """
        sample_query = "INSERT INTO samples (sampleID, studyID) " \
                       "VALUES (%s,%s)"
        self.value_query(sample_query, values)

    def add_meta(self, values):
        """
        Adds rows to the meta table in the PostgreSQL database.

        :param values: List of tuples for meta table
        :return:
        """
        sample_query = "INSERT INTO meta (sampleID, studyID, property, textvalue, numvalue) " \
                       "VALUES (%s,%s,%s,%s,%s)"
        self.value_query(sample_query, values)

    def add_observation(self, values):
        """
        Adds rows to the counts table in the PostgreSQL database.

        :param values: List of tuples for counts table
        :return:
        """
        counts_query = "INSERT INTO counts (studyID,taxon,sampleID,count) " \
                       "VALUES (%s,%s,%s,%s)"
        self.value_query(counts_query, values)
