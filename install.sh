#!/usr/bin/env bash
set -euo pipefail

module_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${module_root}/cmake/transactional-publish.sh"
vcpkg_commit="42f37e46d37baf7b218966434b91116bddf73461"
vcpkg_repository="https://github.com/microsoft/vcpkg.git"
sdk_repository="https://github.com/hweifluids/SCP-SolverSDK.git"
sdk_tag="7ea8900e714c50773d8338d38301fce6c1b073ef"
sdk_version="0.2.0"
architecture="$(uname -m | tr '[:upper:]' '[:lower:]')"
build_dir="${module_root}/.build/linux-${architecture}"
install_dir="${module_root}/release"
vcpkg_installed_dir="${module_root}/.deps/vcpkg_installed"
vcpkg_triplet="x64-linux-release"
clean=0

resolve_component_superproject_root() {
  local component_name="$1"
  local modules_root candidate_root manifest_path component_path expected_mount
  modules_root="$(dirname -- "${module_root}")"
  [[ "$(basename -- "${modules_root}")" == "modules" ]] || return 1
  candidate_root="$(dirname -- "${modules_root}")"
  manifest_path="${modules_root}/components.json"
  [[ -e "${candidate_root}/.git" && -f "${manifest_path}" ]] || return 1
  component_path="$(awk -v wanted="${component_name}" '
    /"name"[[:space:]]*:/ {
      value = $0
      sub(/^[^:]*:[[:space:]]*"/, "", value)
      sub(/".*$/, "", value)
      matched = (value == wanted)
    }
    matched && /"path"[[:space:]]*:/ {
      value = $0
      sub(/^[^:]*:[[:space:]]*"/, "", value)
      sub(/".*$/, "", value)
      print value
      exit
    }
  ' "${manifest_path}")"
  [[ -n "${component_path}" && -d "${modules_root}/${component_path}" ]] || return 1
  expected_mount="$(cd -- "${modules_root}/${component_path}" && pwd -P)"
  [[ "${module_root}" == "${expected_mount}" ]] || return 1
  printf '%s\n' "${candidate_root}"
}

superproject_root="$(resolve_component_superproject_root "SCP-ResamplerIO" || true)"
if [[ -n "${superproject_root}" ]]; then
  vcpkg_root="${superproject_root}/.deps/vcpkg"
else
  vcpkg_root="${module_root}/.deps/vcpkg"
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)
      clean=1
      shift
      ;;
    -h|--help)
      echo "Usage: ./install.sh [--clean]"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

command -v cmake >/dev/null 2>&1 || { echo "CMake 3.24 or newer is required." >&2; exit 1; }
command -v git >/dev/null 2>&1 || { echo "Git is required for the pinned vcpkg bootstrap." >&2; exit 1; }

if [[ ${clean} -eq 1 ]]; then
  build_root="${module_root}/.build"
  [[ "${build_root}" == "${module_root}/.build" ]] || { echo "Refusing to clean outside SCP-ResamplerIO." >&2; exit 1; }
  rm -rf -- "${build_root}"
fi

if [[ ! -x "${vcpkg_root}/vcpkg" || ! -f "${vcpkg_root}/scripts/buildsystems/vcpkg.cmake" ]]; then
  if [[ -e "${vcpkg_root}" ]]; then
    echo "The selected vcpkg root is incomplete: ${vcpkg_root}. Remove or repair it, then retry." >&2
    exit 1
  fi
  mkdir -p "${vcpkg_root}"
  git -C "${vcpkg_root}" init
  git -C "${vcpkg_root}" remote add origin "${vcpkg_repository}"
  git -C "${vcpkg_root}" fetch --depth 1 origin "${vcpkg_commit}"
  git -C "${vcpkg_root}" checkout --detach FETCH_HEAD
  "${vcpkg_root}/bootstrap-vcpkg.sh" -disableMetrics
fi

if [[ ! -d "${vcpkg_root}/.git" ]]; then
  echo "The selected vcpkg root is not a Git checkout: ${vcpkg_root}" >&2
  exit 1
fi
actual_vcpkg_origin="$(git -C "${vcpkg_root}" remote get-url origin 2>/dev/null || true)"
if [[ "${actual_vcpkg_origin}" != "${vcpkg_repository}" ]]; then
  echo "The selected vcpkg checkout must use the pinned origin ${vcpkg_repository}: ${vcpkg_root}" >&2
  exit 1
fi
actual_vcpkg_commit="$(git -C "${vcpkg_root}" rev-parse HEAD 2>/dev/null || true)"
if [[ "${actual_vcpkg_commit,,}" != "${vcpkg_commit}" ]]; then
  echo "The selected vcpkg checkout must be exactly commit ${vcpkg_commit}: ${vcpkg_root}" >&2
  exit 1
fi
vcpkg_status="$(git -C "${vcpkg_root}" status --porcelain 2>/dev/null || true)"
if [[ -n "${vcpkg_status}" ]]; then
  echo "The selected vcpkg checkout must be clean: ${vcpkg_root}" >&2
  exit 1
fi

sdk_sibling_root="${module_root}/../SCP-SolverSDK"
sdk_local_root="${module_root}/.deps/SCP-SolverSDK"
if [[ ! -f "${sdk_sibling_root}/CMakeLists.txt" ]] &&
   [[ ! -f "${sdk_local_root}/CMakeLists.txt" ]] &&
   ! compgen -G "${sdk_sibling_root}/release/*/lib/cmake/SCPSolverSDK/SCPSolverSDKConfig.cmake" >/dev/null &&
   ! compgen -G "${sdk_local_root}/release/*/lib/cmake/SCPSolverSDK/SCPSolverSDKConfig.cmake" >/dev/null; then
  if [[ -e "${sdk_local_root}" ]]; then
    echo "The fixed module-local SDK checkout is incomplete: ${sdk_local_root}. Remove or repair it, then retry." >&2
    exit 1
  fi
  mkdir -p "$(dirname -- "${sdk_local_root}")"
  if ! { git init "${sdk_local_root}" &&
         git -C "${sdk_local_root}" remote add origin "${sdk_repository}" &&
         git -C "${sdk_local_root}" fetch --depth 1 origin "${sdk_tag}" &&
         git -C "${sdk_local_root}" checkout --detach "${sdk_tag}"; }; then
    echo "Could not bootstrap SCP-SolverSDK from ${sdk_repository} revision ${sdk_tag}. The revision must be published before detached installation can succeed." >&2
    exit 1
  fi
fi

sdk_source_root=""
if [[ -f "${sdk_sibling_root}/CMakeLists.txt" ]]; then
  sdk_source_root="${sdk_sibling_root}"
elif [[ -f "${sdk_local_root}/CMakeLists.txt" ]]; then
  sdk_source_root="${sdk_local_root}"
fi
if [[ -n "${sdk_source_root}" ]]; then
  if [[ ! -f "${sdk_source_root}/VERSION" ]] ||
     [[ "$(tr -d '[:space:]' < "${sdk_source_root}/VERSION")" != "${sdk_version}" ]]; then
    echo "SCP-ResamplerIO requires SCP-SolverSDK ${sdk_version} at ${sdk_source_root}." >&2
    exit 1
  fi
fi
if [[ "${sdk_source_root}" == "${sdk_local_root}" ]]; then
  actual_sdk_origin="$(git -C "${sdk_local_root}" remote get-url origin 2>/dev/null || true)"
  if [[ "${actual_sdk_origin}" != "${sdk_repository}" ]]; then
    echo "The module-local SCP-SolverSDK checkout must use origin ${sdk_repository}: ${sdk_local_root}" >&2
    exit 1
  fi
  checked_out_tag="$(git -C "${sdk_local_root}" rev-parse HEAD 2>/dev/null || true)"
  if [[ "${checked_out_tag}" != "${sdk_tag}" ]]; then
    echo "The module-local SCP-SolverSDK checkout must be exactly revision ${sdk_tag}: ${sdk_local_root}" >&2
    exit 1
  fi
fi

configure_args=(
  -S "${module_root}"
  -B "${build_dir}"
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_TRY_COMPILE_CONFIGURATION=Release
  -DCMAKE_INSTALL_PREFIX="${install_dir}"
  -DCMAKE_TOOLCHAIN_FILE="${vcpkg_root}/scripts/buildsystems/vcpkg.cmake"
  -DVCPKG_TARGET_TRIPLET="${vcpkg_triplet}"
  -DVCPKG_HOST_TRIPLET="${vcpkg_triplet}"
  -DVCPKG_OVERLAY_TRIPLETS="${module_root}/cmake/triplets"
  -DVCPKG_INSTALLED_DIR="${vcpkg_installed_dir}"
  -DSCP_RESAMPLER_IO_BUILD_TOOL=ON
)
if command -v ninja >/dev/null 2>&1; then
  configure_args+=( -G Ninja )
fi
cmake "${configure_args[@]}"
cmake --build "${build_dir}" --config Release --parallel

if [[ -z "${install_dir}" || "${install_dir}" == "/" ]]; then
  echo "Refusing to replace an empty or filesystem-root install prefix." >&2
  exit 1
fi

validate_resampler_release_identity() {
  local release_root="$1" identity required tool_path
  local -a tools=()
  if [[ -d "${release_root}/tools" ]]; then
    mapfile -d '' tools < <(find "${release_root}/tools" -maxdepth 1 -type f \
      -name 'Resample_solver_*' -print0)
  fi
  if [[ ${#tools[@]} -ne 1 ]]; then
    echo "Validated release must contain exactly one non-empty dated Resample tool under ${release_root}/tools." >&2
    return 1
  fi
  if [[ ! -s "${tools[0]}" ]]; then
    echo "Validated Resample tool is empty: ${tools[0]}" >&2
    return 1
  fi
  tool_path="${tools[0]}"
  identity="$("${tool_path}" --solver-info)" || return 1
  for required in \
      'streamcenterplus_solver_identity=1' \
      'product=Streamcenter+' \
      'type=tools' \
      'solver=Resample' \
      'version=0.3.0' \
      'cpu=1' \
      'cuda=0'; do
    grep -Fqx -- "${required}" <<<"${identity}" || {
      echo "The staged Resample identity is missing exact line '${required}'." >&2
      return 1
    }
  done
}

stage_parent="${build_dir}/publish-stage"
staged_release="${stage_parent}/release"
[[ "${stage_parent}" == "${build_dir}/publish-stage" ]] || {
  echo "Refusing to stage installation outside the module build tree." >&2
  exit 1
}
rm -rf -- "${stage_parent}"
mkdir -p -- "${stage_parent}"
cleanup_publish_stage() {
  rm -rf -- "${stage_parent}"
}
trap cleanup_publish_stage EXIT

cmake --install "${build_dir}" --config Release --prefix "${staged_release}"
validate_resampler_release_identity "${staged_release}"
validated_identity="$(publish_transactional_directory \
  "${staged_release}" "${install_dir}" validate_resampler_release_identity)"
echo "Validated release content SHA-256: ${validated_identity}"
trap - EXIT
cleanup_publish_stage

tool_path="$(find "${install_dir}/tools" -maxdepth 1 -type f -name 'Resample_solver_*' -print -quit)"

echo "SCP-ResamplerIO installed to ${install_dir}"
echo "Resample tool ready at ${tool_path}"
