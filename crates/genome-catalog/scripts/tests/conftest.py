import sys
from pathlib import Path

# Let the tests import make_genome_list from the scripts directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
