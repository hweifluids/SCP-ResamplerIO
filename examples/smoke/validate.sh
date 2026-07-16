#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
module_root="$(cd -- "${script_directory}/../.." && pwd)"
tool="${module_root}/release/tools/scp_nek5000_smoke"
resolve_component_superproject_root() {
  local component_name="$1" modules_root candidate_root manifest_path component_path expected_mount
  modules_root="$(dirname -- "${module_root}")"
  [[ "$(basename -- "${modules_root}")" == "modules" ]] || return 1
  candidate_root="$(dirname -- "${modules_root}")"
  manifest_path="${modules_root}/components.json"
  [[ -e "${candidate_root}/.git" && -f "${manifest_path}" ]] || return 1
  component_path="$(awk -v wanted="${component_name}" '
    /"name"[[:space:]]*:/ { value=$0; sub(/^[^:]*:[[:space:]]*"/, "", value); sub(/".*$/, "", value); matched=(value==wanted) }
    matched && /"path"[[:space:]]*:/ { value=$0; sub(/^[^:]*:[[:space:]]*"/, "", value); sub(/".*$/, "", value); print value; exit }
  ' "${manifest_path}")"
  [[ -n "${component_path}" && -d "${modules_root}/${component_path}" ]] || return 1
  expected_mount="$(cd -- "${modules_root}/${component_path}" && pwd -P)"
  [[ "${module_root}" == "${expected_mount}" ]] || return 1
  printf '%s\n' "${candidate_root}"
}
superproject_root="$(resolve_component_superproject_root "SCP-ResamplerIO" || true)"
if [[ -n "${superproject_root}" ]]; then
  test_root="${superproject_root}/.tests"
else
  test_root="${module_root}/.tests"
fi
dataset_directory="${test_root}/smoke/SCP-ResamplerIO/data"

if [[ ! -x "${tool}" ]]; then
  echo "The smoke tool is not installed at its fixed location: ${tool}. Run ../../install.sh first." >&2
  exit 2
fi
"${tool}" validate "${dataset_directory}"
