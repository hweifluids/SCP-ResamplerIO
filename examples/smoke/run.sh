#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_root="$(cd -- "${script_directory}/../.." && pwd)"

bash "${module_root}/install.sh"
bash "${script_directory}/prepare_data.sh"
