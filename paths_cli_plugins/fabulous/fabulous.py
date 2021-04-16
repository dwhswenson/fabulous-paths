import os
import click
from fabulous_paths import main
import mdtraj as md
from paths_cli.parameters import INPUT_FILE, MULTI_CV

@click.command(
    "fabulous",
    short_help="Analysis using the FABULOUS framework"
)
@INPUT_FILE.clicked(required=True)
@MULTI_CV.clicked(required=True)   # TODO: option to allow --cv-list file?
@click.option('--ref', type=click.Path(readable=True), required=True,
              help="reference frame for alignment")
@click.option('--keep-atoms', type=str, required=True,
              help="atoms to keep")  # TODO: also allow `.npy` files
@click.option('-c', '--conf', type=click.Path(readable=True), required=True,
              help="FABULOUS configuration file")
@click.option('-n', '--ngen', type=int,
              help="number of generations to run")
@click.option('--results', type=click.Path(writable=True), required=True,
              help="directory to store results in")
@click.option('--label', type=str, required=False, default="1",
              help="label for this run")
def fabulous(input_file, cv, ref, keep_atoms, conf, ngen, results, label):
    storage = INPUT_FILE.get(input_file)
    cvs = MULTI_CV.get(storage, cv)
    ref = md.load(ref)
    main(storage.steps, cvs, ref, keep_atoms, conf, ngen, results, label)


# _mock_command used in docs generation
_mock_command = click.Group(name='openpathsampling',
                            commands={'fabulous': fabulous})


CLI = fabulous
SECTION = "Analysis"
OPS_VERSION = (1, 5)
