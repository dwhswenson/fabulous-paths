import click


@click.command(
    "fabulous",
    short_help="Analysis using the FABULOUS framework"
)
def fabulous():
    print("It works")

CLI = fabulous
SECTION = "Analysis"
OPS_VERSION = (1, 0)
