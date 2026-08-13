#!/usr/bin/env bash
set -euo pipefail

helper="$1"
fixture_root="$2"
source "${helper}"

case "${fixture_root}" in
  */.tests/*/install-publish-linux) ;;
  *)
    echo "Refusing transactional publish test fixture outside .tests/.../install-publish-linux: ${fixture_root}" >&2
    exit 1
    ;;
esac

rm -rf -- "${fixture_root}"
mkdir -p -- "${fixture_root}/release" "${fixture_root}/stage-one"
printf old >"${fixture_root}/release/payload.txt"
printf new-one >"${fixture_root}/stage-one/payload.txt"

validate_one() {
  [[ "$(cat -- "$1/payload.txt")" == "new-one" ]]
}
publish_transactional_directory "${fixture_root}/stage-one" "${fixture_root}/release" validate_one >/dev/null
[[ "$(cat -- "${fixture_root}/release/payload.txt")" == "new-one" ]]

mkdir -p -- "${fixture_root}/stage-mutation"
printf mutation >"${fixture_root}/stage-mutation/payload.txt"
mutated_stage=0
mutate_after_identity() {
  if [[ ${mutated_stage} -eq 0 ]]; then
    printf '%s' '-changed' >>"$1/payload.txt"
    mutated_stage=1
  fi
}
if publish_transactional_directory \
     "${fixture_root}/stage-mutation" "${fixture_root}/release" mutate_after_identity >/dev/null; then
  echo "Staged content mutation was not rejected." >&2
  exit 1
fi
[[ "$(cat -- "${fixture_root}/release/payload.txt")" == "new-one" ]]

mkdir -p -- "${fixture_root}/stage-two"
printf new-two >"${fixture_root}/stage-two/payload.txt"
validation_calls=0
fail_after_swap() {
  validation_calls=$((validation_calls + 1))
  [[ "$(cat -- "$1/payload.txt")" == "new-two" ]]
  [[ ${validation_calls} -lt 3 ]]
}
if publish_transactional_directory "${fixture_root}/stage-two" "${fixture_root}/release" fail_after_swap >/dev/null; then
  echo "Injected post-swap failure was not reported." >&2
  exit 1
fi
[[ "$(cat -- "${fixture_root}/release/payload.txt")" == "new-one" ]]
echo "Transactional Linux publish and rollback test passed."
