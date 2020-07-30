# _masq_

Microbial association network SQL API.

This API contains drivers for interacting with a PostgreSQL database.
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


