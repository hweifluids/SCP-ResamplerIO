#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cmake "-DSMOKE_ROOT=${script_directory}" -P "${script_directory}/prepare-data.cmake"
bash "${script_directory}/validate.sh"
