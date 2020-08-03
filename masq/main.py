"""
This file contains parsers and functions that call on other functionality defined
in the rest of masq's scripts directory.
The command line interface is intended to be called sequentially;
files are written to disk as intermediates,
while a settings file is used to transfer logs and other information
between the modules. These modules are contained in this file.
This modular design allows users to leave out parts of masq that are not required,
and reduces the number of parameters that need to be defined in the function calls.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import sys
import argparse
from ast import literal_eval
from pbr.version import VersionInfo
from masq.scripts.sq4biom import import_biom
from masq.scripts.io import import_networks
from masq.scripts.utils import setup_database
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


def masq(masq_args):
    """
    Main function for running masq.
    Accepts a dictionary of arguments from the argument parser
    and calls the appropriate module function.

    :param masq_args: Arguments.
    :return:
    """
    if masq_args['version']:
        info = VersionInfo('anuran')
        logger.info('Version ' + info.version_string())
        sys.exit(0)
    # unpack args
    config = masq_args['config']
    database = masq_args['config']
    networks = masq_args['networks']
    bioms = masq_args['bioms']
    username = masq_args['username']
    password = masq_args['password']
    host = masq_args['host']
    mapping = masq_args['mapping']
    if masq_args['mapping']:
        try:
            with open(masq_args['mapping'], 'r') as file:
                contents = file.read()
                mapping = literal_eval(contents)
        except (ValueError, TypeError):
            logger.warning("Mapping file could not be imported,\n"
                           "and will be ignored. ")
    if masq_args['create']:
        logger.info('Setting up tables in PostgreSQL database. ')
        setup_database(config=config, create=True,
                       host=host, database=database,
                       username=username, password=password)
    elif masq_args['delete']:
        setup_database(config=config, create=False,
                       host=host, database=database,
                       username=username, password=password)
        logger.info("Deleted tables in PostgreSQL database.")
    if masq_args['bioms']:
        logger.info('Importing BIOM files... ')
        for biom in bioms:
            import_biom(location=biom, mapping=mapping,
                        config=config,
                        host=host, database=database,
                        username=username, password=password)
    if masq_args['networks']:
        logger.info('Importing network files...')
        for network in networks:
            import_networks(location=network, mapping=mapping,
                            config=config,
                            host=host, database=database,
                            username=username, password=password)
    logger.info('Completed tasks! ')


masq_parser = argparse.ArgumentParser(description='masq API')
masq_parser.add_argument('-c', '--config',
                         dest='config',
                         help='Config file containing '
                              'default settings. ',
                         default='database.ini',
                         type=str)
masq_parser.add_argument('-d', '--database',
                         dest='database',
                         help='Name of PostgreSQL database. '
                              'Can also be specified in config file.',
                         default=None,
                         type=str)
masq_parser.add_argument('-ht', '--host',
                         dest='host',
                         help='Host location for PostgreSQL database. '
                              'Can also be specified in config file.',
                         default=None,
                         type=str)
masq_parser.add_argument('-u', '--username',
                         dest='username',
                         help='Username for PostgreSQL database. '
                              'Can also be specified in config file.',
                         default=None,
                         type=str)
masq_parser.add_argument('-p', '--password',
                         dest='password',
                         help='Password for PostgreSQL database. '
                              'Can also be specified in config file.',
                         default=None,
                         type=str)
masq_parser.add_argument('-p', '--password',
                         dest='password',
                         help='Password for PostgreSQL database. '
                              'Can also be specified in config file.',
                         default=None,
                         type=str)
masq_parser.add_argument('-c', '--create',
                         dest='create',
                         help='If flagged, sets up tables in PostgreSQL database. ',
                         default=False,
                         action='store_true')
masq_parser.add_argument('-del', '--delete',
                         dest='delete',
                         help='If flagged, deletes tables in PostgreSQL database. ',
                         default=False,
                         action='store_true')
masq_parser.add_argument('-b', '--bioms',
                         dest='bioms',
                         help='One or more filenames or folders containing BIOM files.\n'
                              'These are imported in the PostgreSQL database. ',
                         default=None,
                         type=list)
masq_parser.add_argument('-n', '--networks',
                         dest='networks',
                         help='One or more filenames or folders containing network files.\n'
                              'These are imported in the PostgreSQL database. ',
                         default=None,
                         type=list)
masq_parser.add_argument('-version', '--version',
                         dest='version',
                         required=False,
                         help='Version number.',
                         action='store_true',
                         default=False)


def main():
    options = masq_parser.parse_args()
    masq(vars(options))


if __name__ == '__main__':
    main()
