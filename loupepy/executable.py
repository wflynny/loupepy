import platform
import pkg_resources
from hashlib import md5
from pathlib import Path
from tempfile import mkstemp
from urllib.request import urlretrieve


BINARY_ARTIFACTS = dict(
    Linux = dict(
        url = "https://github.com/10XGenomics/loupeR/releases/download/v1.1.0/louper-linux-x64",
        md5 = "99903df7a3bc7b1b06d7509ddddf9a13",
        target = "louper"
    ),
    Darwin = dict(
        url = "https://github.com/10XGenomics/loupeR/releases/download/v1.1.0/louper-macos-x64",
        md5 = "bf4ff652b88e0b9a88fb306b11a9c066",
        target = "louper"
    ),
    Windows = dict(
        url = "https://github.com/10XGenomics/loupeR/releases/download/v1.1.0/louper-windows-x64.exe",
        md5 = "f40833260e3d4c14d8534a1f3349096d",
        target = "louper.exe"
    )
)


def setup_executable(executable_path: Path = None):
    if executable_path and executable_path.exists():
        return executable_path

    if find_executable() is None:
        print("Downloading executable...")
        executable_path = download_artifact()
        return executable_path

    

def find_executable():
    default = default_executable_path()
    if default:
        return default
    bundled = artifact_is_bundled()
    if bundled:
        return bundled


def get_artifact():
    os_name = platform.system()
    os_dict = BINARY_ARTIFACTS.get(os_name, None)
    if os_dict is None:
        raise Exception(f"Invalid OS {os_name}")
    return os_dict


def artifact_is_bundled():
    artifact_dict = get_artifact()
    target = artifact_dict["target"]
    putative_exe = Path(pkg_resources.resource_string("loupepy", f"bin/{target}"))
    return putative_exe if putative_exe.exists() else None


def get_data_dir(data_dir=None) -> Path:
    if data_dir is None:
        data_dir = Path.home() / ".config" / "loupepy"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def default_executable_path():
    os_dict = get_artifact()
    return get_data_dir() / os_dict["target"]


def download_artifact() -> None:
    os_dict = get_artifact()
    
    tmpfile = Path(mkstemp()[1])
    urlretrieve(os_dict["url"], tmpfile)

    dest = get_data_dir() / os_dict["target"]
    tmpfile.rename(dest)
    dest.chmod(0o0755)

    if not dest.exists():
        raise FileNotFoundError(f"Download and move failed for path {dest}")

    hasher = md5()
    with open(dest, "rb") as b:
        hasher.update(b.read())
        digest = hasher.hexdigest()
    if digest != os_dict["md5"]:
        raise Exception(f"Executable {dest} md5sum doesn't match expected\nfound: {digest}\ntarget: {os_dict['md5']}")

    return dest


