"""
Database download module for PlasAnn
Downloads databases from Zenodo with Google Drive fallback
"""

import requests
import zipfile
import hashlib
import gdown
import pandas as pd
import subprocess
from pathlib import Path

# ⚠️ REPLACE "XXXXXXX" WITH YOUR ACTUAL ZENODO RECORD ID ⚠️
ZENODO_CONFIG = {
    "record_id": "15583460",  # ← Your actual Zenodo record ID
    "version": "v2",
    "zip_filename": "plasann_databases_v1.1.0.zip"
}

# Fallback source: Google Drive  
GDRIVE_CONFIG = {
    "folder_id": "14jAiNrnsD7p0Kje--nB23fq_zTTE9noz"
}

def get_zenodo_record_info(record_id):
    """Get information about a Zenodo record using the API"""
    try:
        url = f"https://zenodo.org/api/records/{record_id}"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"❌ Could not fetch Zenodo record info: {e}")
        return None

def download_from_zenodo():
    """Download databases from Zenodo (primary method)"""
    print("🏛️ Downloading from Zenodo (primary source)...")
    
    record_id = ZENODO_CONFIG["record_id"]
    zip_filename = ZENODO_CONFIG["zip_filename"]
    
    if record_id == "XXXXXXX":
        raise Exception("❌ Zenodo record ID not configured! Please update ZENODO_CONFIG['record_id']")
    
    # Get record info via API
    print(f"🔍 Fetching Zenodo record {record_id}...")
    record_data = get_zenodo_record_info(record_id)
    
    if not record_data:
        raise Exception("Could not fetch Zenodo record information")
    
    # Find the database file
    target_file = None
    for file_info in record_data.get("files", []):
        if file_info["key"] == zip_filename:
            target_file = file_info
            break
    
    if not target_file:
        available_files = [f["key"] for f in record_data.get("files", [])]
        raise Exception(f"Database file '{zip_filename}' not found. Available files: {available_files}")
    
    download_url = target_file["links"]["self"]
    file_size = target_file["size"]
    checksum = target_file.get("checksum", "")
    
    # Download with progress
    print(f"⬇️ Downloading {zip_filename} ({file_size // 1024 // 1024}MB)")
    print(f"   URL: {download_url}")
    
    response = requests.get(download_url, stream=True)
    response.raise_for_status()
    
    zip_path = Path("databases.zip")
    downloaded = 0
    
    with open(zip_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                if file_size > 0:
                    percent = (downloaded / file_size) * 100
                    print(f"\r   Progress: {percent:.1f}% ({downloaded//1024//1024}MB/{file_size//1024//1024}MB)", end='')
    
    print(f"\n✅ Zenodo download completed")
    
    # Verify checksum if available
    if checksum and checksum.startswith("md5:"):
        expected_hash = checksum.split(":", 1)[1]
        print("🔍 Verifying file integrity...")
        if verify_md5_checksum(zip_path, expected_hash):
            print("✅ File integrity verified")
        else:
            zip_path.unlink()
            raise Exception("❌ File checksum verification failed!")
    
    return zip_path

def verify_md5_checksum(file_path, expected_hash):
    """Verify MD5 checksum of downloaded file"""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest() == expected_hash

def download_from_gdrive():
    """Download databases from Google Drive (fallback)"""
    print("📁 Downloading from Google Drive (fallback)...")
    
    db_folder = Path("Database")
    db_folder.mkdir(exist_ok=True)
    
    folder_id = GDRIVE_CONFIG["folder_id"]
    
    # Try multiple times with gdown (it can be flaky)
    max_retries = 3
    for attempt in range(max_retries):
        try:
            print(f"   Attempt {attempt + 1}/{max_retries}...")
            gdown.download_folder(
                f"https://drive.google.com/drive/folders/{folder_id}",
                output=str(db_folder),
                quiet=False,
                use_cookies=False,
                remaining_ok=True
            )
            print("✅ Google Drive download completed")
            return None  # Direct extraction, no ZIP
            
        except Exception as e:
            print(f"   ❌ Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                print("   🔄 Retrying...")
            else:
                raise

def extract_databases(zip_path):
    """Extract database ZIP file"""
    if zip_path is None:
        return  # Already extracted (Google Drive method)
    
    print("📦 Extracting database files...")
    db_folder = Path("Database")
    db_folder.mkdir(exist_ok=True)
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for member in zip_ref.namelist():
            # Skip directory entries and hidden files
            if member.endswith('/') or member.startswith('.'):
                continue
                
            # Handle different ZIP structures
            if member.startswith('Database/'):
                # Remove 'Database/' prefix if present in ZIP
                target_filename = member[9:]  
            else:
                target_filename = member
            
            if target_filename:  # Skip empty filenames
                target_path = db_folder / target_filename
                target_path.parent.mkdir(parents=True, exist_ok=True)
                
                with zip_ref.open(member) as source, open(target_path, 'wb') as target:
                    target.write(source.read())
    
    # Clean up ZIP
    zip_path.unlink()
    print("✅ Database extraction completed")

def verify_database():
    """Verify all essential database files are present"""
    db_folder = Path("Database")
    essential_files = [
        "Database.csv",
        "orit.fna", 
        "oriT_TNcentral.fasta",
        "plasmidfinder.fasta",
        "tncentral_cleaned.fa",
        "transposon.fasta", 
        "Rfam.cm",
        "oric.fna",
        "uniprot_plasmids.tsv"  # Added UniProt TSV file
    ]
    
    missing_files = []
    for filename in essential_files:
        file_path = db_folder / filename
        if not file_path.exists():
            missing_files.append(filename)
    
    if missing_files:
        raise Exception(f"Missing essential database files: {missing_files}")
    
    print(f"✅ All {len(essential_files)} essential database files verified")
    
    # Show file sizes for confirmation
    total_size = 0
    for filename in essential_files:
        file_path = db_folder / filename
        if file_path.exists():
            size = file_path.stat().st_size
            total_size += size
    
    print(f"   📊 Total database size: {total_size // 1024 // 1024}MB")

def download_database():
    """Main function to download databases with fallback"""
    db_folder = Path("Database")
    db_folder.mkdir(exist_ok=True)
    
    version_file = db_folder / ".version"
    current_version = ZENODO_CONFIG["version"]
    
    # First, check if all essential files are present
    try:
        verify_database()  # Check if all files exist
        
        # If files exist but no version file, assume it's current version
        if not version_file.exists():
            print("✅ Database files detected, marking as current version...")
            with open(version_file, 'w') as f:
                f.write(current_version)
            print("✅ Database already present and up-to-date.")
            return
        
        # If version file exists, check version
        with open(version_file, 'r') as f:
            installed_version = f.read().strip()
        if installed_version == current_version:
            print("✅ Database already present and up-to-date.")
            return
        else:
            print(f"⚠️ Database version mismatch: {installed_version} → {current_version}")
            print("   Updating to latest version...")
            
    except Exception as e:
        print(f"⚠️ Database missing or incomplete: {e}")
        print("   Downloading fresh database...")
    
    print(f"⬇️ Downloading PlasAnn databases (version {current_version})")
    
    # Try Zenodo first, then Google Drive
    download_sources = [
        ("Zenodo", download_from_zenodo),
        ("Google Drive", download_from_gdrive)
    ]
    
    for source_name, download_func in download_sources:
        try:
            print(f"\n🔄 Trying {source_name}...")
            zip_path = download_func()
            extract_databases(zip_path)
            verify_database()
            
            # Write version file
            with open(version_file, 'w') as f:
                f.write(current_version)
            
            print(f"✅ Successfully downloaded databases from {source_name}")
            print(f"   📁 Database location: {db_folder.absolute()}")
            return
            
        except Exception as e:
            print(f"❌ {source_name} failed: {e}")
            # Clean up any partial downloads
            if zip_path and zip_path.exists():
                zip_path.unlink()
            continue
    
    raise RuntimeError("❌ Failed to download databases from all sources!")

def prepare_blast_database():
    """Prepare BLAST database from downloaded Database.csv"""
    db_csv = Path("Database/Database.csv")
    
    if not db_csv.exists():
        print("❌ Database.csv not found. Running download_database()...")
        download_database()
        
        # Check again after download
        if not db_csv.exists():
            raise FileNotFoundError("Database.csv still not found after download!")

    blast_folder = Path("database_blast")
    blast_folder.mkdir(exist_ok=True)
    fasta_path = blast_folder / "translations.fasta"

    print("🧬 Converting Database.csv to FASTA format...")
    try:
        df = pd.read_csv(db_csv)
        with open(fasta_path, "w") as fasta_out:
            for idx, row in df.iterrows():
                sequence = row.get("Translation", "")
                if pd.notnull(sequence) and sequence.strip():
                    fasta_out.write(f">{idx}\n{sequence.strip()}\n")
        print(f"✅ Database FASTA written with {len(df)} sequences.")
    except Exception as e:
        raise Exception(f"Failed to process Database.csv: {e}")

    # Build BLAST database
    db_prefix = blast_folder / "translations_db"
    if not (blast_folder / "translations_db.pin").exists():
        print("🛠️ Building BLAST database...")
        try:
            subprocess.run([
                "makeblastdb",
                "-in", str(fasta_path),
                "-dbtype", "prot",
                "-out", str(db_prefix)
            ], check=True, capture_output=True, text=True)
            print("✅ BLAST database created.")
        except subprocess.CalledProcessError as e:
            raise Exception(f"Failed to create BLAST database: {e.stderr}")
        except FileNotFoundError:
            raise Exception("makeblastdb not found! Please install BLAST+ tools.")
    else:
        print("✅ BLAST database already exists.")

def get_database_info():
    """Get information about the current database installation"""
    db_folder = Path("Database")
    version_file = db_folder / ".version"
    
    if not version_file.exists():
        return "No database installed"
    
    try:
        with open(version_file, 'r') as f:
            version = f.read().strip()
        
        # Count files and calculate size
        db_files = list(db_folder.glob("*"))
        total_size = sum(f.stat().st_size for f in db_files if f.is_file())
        
        return {
            "version": version,
            "files": len(db_files),
            "size_mb": total_size // 1024 // 1024,
            "zenodo_record": ZENODO_CONFIG["record_id"],
            "location": str(db_folder.absolute())
        }
    except Exception as e:
        return f"Error reading database info: {e}"

def force_redownload():
    """Force redownload of databases (removes version file)"""
    db_folder = Path("Database")
    version_file = db_folder / ".version"
    
    if version_file.exists():
        version_file.unlink()
        print("🔄 Forcing database redownload...")
    
    download_database()

# For testing
if __name__ == "__main__":
    print("🧪 Testing PlasAnn database download...")
    try:
        download_database()
        print("✅ Database download test successful!")
        
        info = get_database_info()
        print(f"📊 Database info: {info}")
        
    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()