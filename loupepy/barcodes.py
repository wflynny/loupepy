import pkg_resources
from pathlib import Path
import pandas as pd


def load_valid_barcodes():
    barcode_file = Path(pkg_resources.resource_filename("loupepy", "barcodes/3M-3pgex-may-2023.txt.gz"))
    bcs = pd.read_table(barcode_file, index_col=0, header=None).index
    bcs.name = None
    return bcs

VALID_BARCODES = load_valid_barcodes()
