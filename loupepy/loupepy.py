from pathlib import Path
from subprocess import run
from tempfile import mkstemp
from anndata import AnnData

from .hdf5 import create_hdf5
from .executable import find_executable, setup_executable


def create_loupe_from_anndata(
    adata,
    output_dir = None,
    output_name = None,
    executable_path = None,
    overwrite = False
):

    executable_path = setup_executable(executable_path)
    if not isinstance(adata, AnnData):
       raise Exception("Input object is not an AnnData object")

    h5file = Path(mkstemp[1])
    create_hdf5(adata, h5file)
    create_loupe(
        h5file, 
        output_dir=output_dir,
        output_name=output_name,
        executable_path=executable_path,
        force=overwrite
    )
    h5file.unlink()
    

def create_loupe(
    h5path: Path,
    output_dir = None,
    output_name:str = None, 
    executable_path:Path = None,
    force = False
):
    if not output_name:
        output_name = "converted"
    
    if not output_dir:
        output_dir = "."
    loupe_path = Path(f"{output_dir}/{output_name}.cloupe").resolve()

    if loupe_path.exists() and (not force):
        raise Exception(f"Loupe file '{loupe_path}' already exists. Set `force=True` to overwrite")
    
    args = [
        "create", 
        f"--input={h5path}", 
        f"--output={loupe_path}"
        "--force" if force else ""
    ]

    if not executable_path:
        executable_path = find_executable()
        if not executable_path:
            raise FileNotFoundError("Could not find a valid louper executable. Run `setup_executable()`.")
    if not executable_path.exists():
        raise FileNotFoundError(f"Executable {executable_path} does not exist.")

    print(f"Running command: {executable_path} {' '.join(args)}")

    res = run([str(executable_path)] + args, shell=True)
    if res.returncode > 0:
        raise Exception(f"Louper executable failed: status code {res.returncode}")