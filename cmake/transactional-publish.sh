#!/usr/bin/env bash

transactional_directory_identity() {
  local root="$1" relative hash size
  [[ -d "${root}" ]] || { echo "Publish staging directory does not exist: ${root}" >&2; return 1; }
  [[ ! -L "${root}" ]] || { echo "Publish staging root must not be a symbolic link: ${root}" >&2; return 1; }
  if find "${root}" -type l -print -quit | grep -q .; then
    echo "Publish staging content must not contain symbolic links: ${root}" >&2
    return 1
  fi
  local identity_input identity_parent
  identity_parent="$(dirname -- "${root}")"
  identity_input="$(mktemp "${identity_parent}/.transactional-identity.XXXXXX")" || return 1
  : >"${identity_input}"
  local count=0
  while IFS= read -r -d '' file; do
    relative="${file#"${root}"/}"
    hash="$(sha256sum -- "${file}" | awk '{print $1}')" || { rm -f -- "${identity_input}"; return 1; }
    size="$(wc -c <"${file}" | tr -d '[:space:]')"
    printf '%s:%s|%s|%s\n' "${#relative}" "${relative}" "${size}" "${hash}" >>"${identity_input}"
    count=$((count + 1))
  done < <(find "${root}" -type f -print0 | sort -z)
  if [[ ${count} -eq 0 ]]; then
    rm -f -- "${identity_input}"
    echo "Publish staging directory contains no files: ${root}" >&2
    return 1
  fi
  sha256sum -- "${identity_input}" | awk '{print $1}'
  local status=${PIPESTATUS[0]}
  rm -f -- "${identity_input}"
  return "${status}"
}

publish_transactional_directory() {
  local stage="$1" destination="$2" validator="$3"
  [[ -d "${stage}" ]] || { echo "Publish staging directory does not exist: ${stage}" >&2; return 1; }
  local parent leaf token candidate backup failed identity moved_old=0 stage_path destination_path
  parent="$(cd -- "$(dirname -- "${destination}")" && pwd -P)" || return 1
  leaf="$(basename -- "${destination}")"
  stage_path="$(cd -- "${stage}" && pwd -P)" || return 1
  destination_path="${parent}/${leaf}"
  [[ -n "${leaf}" && "${leaf}" != "." && "${leaf}" != ".." ]] || {
    echo "Unsafe publish destination: ${destination}" >&2
    return 1
  }
  case "${stage_path}/" in
    "${destination_path}/"|"${destination_path}/"*)
      echo "Publish staging and destination directories are not safely separated." >&2
      return 1
      ;;
  esac
  case "${destination_path}/" in
    "${stage_path}/"*)
      echo "Publish staging and destination directories are not safely separated." >&2
      return 1
      ;;
  esac
  if [[ -e "${destination}" && ( ! -d "${destination}" || -L "${destination}" ) ]]; then
    echo "Existing publish destination must be an ordinary directory: ${destination}" >&2
    return 1
  fi
  token="$$.$RANDOM.$RANDOM"
  candidate="${parent}/.${leaf}.publish.${token}"
  backup="${parent}/.${leaf}.rollback.${token}"
  failed="${parent}/.${leaf}.failed.${token}"
  [[ ! -e "${candidate}" && ! -e "${backup}" && ! -e "${failed}" ]] || {
    echo "Transactional publish path already exists." >&2
    return 1
  }

  identity="$(transactional_directory_identity "${stage}")" || return 1
  "${validator}" "${stage}" || return 1
  mv -- "${stage}" "${candidate}" || return 1
  if [[ "$(transactional_directory_identity "${candidate}")" != "${identity}" ]] ||
     ! "${validator}" "${candidate}"; then
    rm -rf -- "${candidate}"
    return 1
  fi

  if [[ -e "${destination}" ]]; then
    mv -- "${destination}" "${backup}" || { rm -rf -- "${candidate}"; return 1; }
    moved_old=1
  fi
  if ! mv -- "${candidate}" "${destination}" ||
     [[ "$(transactional_directory_identity "${destination}")" != "${identity}" ]] ||
     ! "${validator}" "${destination}"; then
    if [[ -e "${destination}" ]]; then
      mv -- "${destination}" "${failed}" || true
    fi
    if [[ ${moved_old} -eq 1 && -e "${backup}" ]]; then
      mv -- "${backup}" "${destination}" || true
    fi
    rm -rf -- "${candidate}" "${failed}"
    return 1
  fi
  rm -rf -- "${backup}" "${failed}"
  printf '%s\n' "${identity}"
}
