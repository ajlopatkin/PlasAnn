"""
Command-line interface for PlasAnn
"""

from .annotate_plasmid import main

def run_main():
    """
    Entry point for the CLI command.
    This function is a simple wrapper around the actual main function.
    """
    main()