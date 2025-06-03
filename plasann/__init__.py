__version__ = "1.1.1"

from .annotate import main as annotate_main, cli_main
# Update these imports to match the new function names
from .annotate import run_on_single_fasta, run_on_single_genbank

def main():
    """Entry point for PlasAnn command"""
    cli_main()

# Keep the old names for backward compatibility
run_on_fasta = run_on_single_fasta
run_on_genbank = run_on_single_genbank