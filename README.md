`loupepy` creates a 10x Genomics Loupe file from a Scanpy AnnData object. This package aims to be the Python equivalent of [loupeR](https://github.com/10XGenomics/loupeR), and as such, *only single-cell gene expression datasets are supported.*

## How to Use
Use as simply as
```{python}
from loupepy import create_loupe_from_anndata

create_loupe_from_anndata(adata)
```

## Depdencies
This package requires
```
- anndata
- pandas
- numpy
- scipy
```