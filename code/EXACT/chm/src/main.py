import argparse
import sys
from libchm import chm

"""
Rational arithmetic implementation of the Convex Hull algorithm for enumeration of vertices in the projection of a polytope onto the dimensions of interest. 

@author: Sarah Noel Galleguillos, Norbert Auer
@email_address: sarah.galleguillos@boku.ac.at, norbert.auer@boku.ac.at
"""


def start(args):
    """
    Prints or save view
    :param args: argparse
    :return: None
    """
    chm.compute_CH(args.name, args.dimensions, args.output)


# Create Argumentparser
parser = argparse.ArgumentParser(prog="chm", description="Returns the vertices of the convex hull for the projection of a polytope onto the dimensions of interest")

# Add an optinal argument to the main parser
parser.add_argument("name", help="File name without extension, must be the same for all three files: filename.d with Rmin Rmax for each variable, filename.r with variable names and filename.S with integer Stoichiometric matrix", type=str)
parser.add_argument("-d", "--dimensions", help="indices of the dimensions onto which the projection should be computed in Python zero based numbering", type=int, nargs="+", required=True)
parser.add_argument("-o", "--output", help="List of Extreme Points of the projection", type=argparse.FileType('w'), default=sys.stdout)

# This function is called automatically if no subparser is used
parser.set_defaults(func=start)


# Main Function
def main():
    args = parser.parse_args()

    # Run default function set by argparse choice (e.g. func=view)
    args.func(args)


if __name__ == "__main__":
    # Runs the main function if called directly
    main()
