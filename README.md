# _masq_

Microbial association network SQL API.

This API contains drivers for interacting with a [PostgreSQL](https://www.postgresql.org/) database.
Many of these drivers have functions that allow for BIOM files,
network files and other microbiome-related files to be ported to the database.

Contact the author at lisa.rottjers (at) kuleuven.be. Your feedback is much appreciated!
This version is still in early beta and has been tested for Python 3.6.

## Getting Started

First set up a [virtual environment](https://docs.python-guide.org/dev/virtualenvs/) and make sure it uses Python 3:
```
virtualenv venv
# Linux
source venv/bin/activate

# Windows
venv/Scripts/activate

# Once you are done with masq:
deactivate
```

To install _masq_, run:
```
pip install git+https://github.com/ramellose/masq.git
```

If you have Python 3 and 2 installed and run Windows, you may need this command:
```
python3 -m pip install git+https://github.com/ramellose/masq.git
```

At the moment, the CLI is not very well documented, since _masq_ is in early development.
You can run the _masq_ script and read the help docs with the following command.

```
masq -h
```

_masq_ requires a running instance of PostgreSQL to connect to.
It will first set up a database structure in line with the database model used by [_mako_](https://github.com/ramellose/mako).
Then it can upload BIOM files and network files to the PostgreSQL database.
Additionally, the _metastats_ and _netstats_ modules contain scripts with queries that can carry out some operations on the PostgreSQL database.
These are currently not supported by the CLI, only by the API.

For viewing PostgreSQL databases, I recommend [HeidiSQL](https://www.heidisql.com/).

If you have both Python 2.7 and Python 3 installed, you may need to change the command to this:
```
python3 -m pip install git+https://github.com/ramellose/masq.git
```

### Contributions

This software is still in early alpha. Any feedback or bug reports will be much appreciated!

## Authors

* **Lisa Röttjers** - [ramellose](https://github.com/ramellose)

See also the list of [contributors](https://github.com/ramellose/manta/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.txt](LICENSE.txt) file for details


