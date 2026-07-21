#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version="5.78-109.0"
base_url="https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/$version"
archive_name="interproscan-${version}-64-bit.tar.gz"
download_dir="$project_dir/tools/downloads"
install_dir="$project_dir/tools/interproscan-$version"

mkdir -p "$download_dir" "$project_dir/tools"

if [[ -x "$install_dir/interproscan.sh" ]]; then
    echo "InterProScan already installed: $install_dir"
    exit 0
fi

curl \
    --fail \
    --location \
    --continue-at - \
    --retry 5 \
    --retry-delay 5 \
    --output "$download_dir/$archive_name" \
    "$base_url/$archive_name"

curl \
    --fail \
    --location \
    --retry 5 \
    --output "$download_dir/$archive_name.md5" \
    "$base_url/$archive_name.md5"

(
    cd "$download_dir"
    md5sum --check "$archive_name.md5"
)

tar -xzf "$download_dir/$archive_name" -C "$project_dir/tools"

if [[ ! -x "$install_dir/interproscan.sh" ]]; then
    echo "InterProScan executable missing after extraction: $install_dir/interproscan.sh" >&2
    exit 1
fi

echo "Installed InterProScan $version at $install_dir"
