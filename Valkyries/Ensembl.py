import pyensembl
import subprocess
import Valkyries.Utilities as Utilities

__author__ = "Dennis A. Simpson"
__version__ = "0.1.2"
__release__ = 93


def initialize(args, log):
    """
    This will initialize pyensembl to the database version below when called.
    :param args:
    :param log:
    :return:
    """

    species = args.Species
    if args.Species == "Mouse" or args.Species == "mus_musculus":
        species = "mus_musculus"
    elif args.Species == "Human" or args.Species == "homo_sapiens":
        species = "Homo_sapiens"
    else:
        log.warning("Unknown species.  Attempting to use {} anyway.".format(args.Species))

    log.info("Initializing Pyensembl --release {} --species {}.  This may take a while.".format(__release__, species))

    t = subprocess.run("pyensembl install --release {} --species {}".format(__release__, species),
                       shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    if t.returncode:
        for line in Utilities.byte_array_to_string(t.stderr).split("\n"):
            log.error(line)
        raise SystemExit(1)
    else:
        for line in Utilities.byte_array_to_string(t.stdout).split("\n"):
            log.debug(line)

    log.info("Pyensembl Initialized")

    try:
        log.debug("Creating Ensembl Data Object.")
        ensembl_data = pyensembl.EnsemblRelease(species=species)
    except ValueError:
        log.error("Pyensembl not initialized correctly")
        raise SystemExit(1)

    return ensembl_data
