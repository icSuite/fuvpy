"""Download the fuvpy validation data from Zenodo.

The archive contains complete IMAGE/FUV orbit sequences and historical
NetCDF background products used by ``scripts/compare_models.py``. It is kept
outside Git because the compressed archive is approximately 5 GB.
"""

import hashlib
import json
from pathlib import Path
import shutil
import urllib.request
import zipfile


REPOSITORY = Path(__file__).resolve().parents[1]
DATA_DIRECTORY = REPOSITORY / "data"
ARCHIVE = REPOSITORY / "data.zip"
EXTRACT_DIRECTORY = REPOSITORY / ".data_extract"
ZENODO_RECORD = "https://zenodo.org/api/records/22828032"


def get_archive_information():
    """Return the name, URL, size, and checksum of the record's ZIP file."""

    with urllib.request.urlopen(ZENODO_RECORD, timeout=60) as response:
        record = json.load(response)

    archives = [file for file in record["files"] if file["key"].lower().endswith(".zip")]
    if len(archives) != 1:
        raise RuntimeError(f"Expected one ZIP archive in Zenodo record 22828032, found {len(archives)}")

    archive = archives[0]
    algorithm, checksum = archive["checksum"].split(":", maxsplit=1)
    if algorithm.lower() != "md5":
        raise RuntimeError(f"Unsupported Zenodo checksum algorithm: {algorithm}")

    return archive["key"], archive["links"]["self"], int(archive["size"]), checksum


def download_archive(url, expected_size, expected_checksum):
    """Download the archive while reporting progress and calculating MD5."""

    digest = hashlib.md5()
    downloaded = 0

    with urllib.request.urlopen(url, timeout=60) as response, ARCHIVE.open("wb") as output:
        while True:
            block = response.read(8 * 1024 * 1024)
            if not block:
                break

            output.write(block)
            digest.update(block)
            downloaded += len(block)
            print(f"\rDownloaded {downloaded / 1e9:.2f} of {expected_size / 1e9:.2f} GB", end="", flush=True)

    print()
    if downloaded != expected_size:
        raise RuntimeError(f"Incomplete download: expected {expected_size} bytes, received {downloaded}")
    if digest.hexdigest() != expected_checksum:
        raise RuntimeError("Downloaded archive failed its Zenodo MD5 checksum")


def extract_archive():
    """Safely extract the archive and move its data folder into the repository."""

    if EXTRACT_DIRECTORY.exists():
        shutil.rmtree(EXTRACT_DIRECTORY)
    EXTRACT_DIRECTORY.mkdir()

    try:
        with zipfile.ZipFile(ARCHIVE) as archive:
            extraction_root = EXTRACT_DIRECTORY.resolve()
            for member in archive.infolist():
                destination = (EXTRACT_DIRECTORY / member.filename).resolve()
                if extraction_root not in destination.parents and destination != extraction_root:
                    raise RuntimeError(f"Unsafe path in ZIP archive: {member.filename}")
            archive.extractall(EXTRACT_DIRECTORY)

        archived_data = EXTRACT_DIRECTORY / "data"
        if archived_data.is_dir():
            shutil.move(str(archived_data), str(DATA_DIRECTORY))
        elif (EXTRACT_DIRECTORY / "idl_files").is_dir() and (EXTRACT_DIRECTORY / "legacy_nc_files").is_dir():
            EXTRACT_DIRECTORY.rename(DATA_DIRECTORY)
        else:
            raise RuntimeError("The downloaded archive does not contain the expected data directory")
    finally:
        if EXTRACT_DIRECTORY.exists():
            shutil.rmtree(EXTRACT_DIRECTORY)


def main():
    if DATA_DIRECTORY.exists():
        raise FileExistsError(f"{DATA_DIRECTORY} already exists; move or remove it before downloading the validation data")

    archive_name, url, size, checksum = get_archive_information()
    print(f"Downloading {archive_name} from Zenodo record 22828032")

    try:
        download_archive(url, size, checksum)
        print("Extracting archive...")
        extract_archive()
    finally:
        if ARCHIVE.exists():
            ARCHIVE.unlink()

    print(f"Validation data installed in {DATA_DIRECTORY}")


if __name__ == "__main__":
    main()
