#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
dest="${1:-$here/vendor/NGA2}"
if [[ -e "$dest" ]]; then echo "Refusing existing target: $dest" >&2; exit 2; fi
git clone --no-checkout https://github.com/desjardi/NGA2.git "$dest"
git -C "$dest" checkout --detach 28a58d937af89dd506456048be1cad58f6d383a9
test "$(git -C "$dest" rev-parse HEAD)" = "28a58d937af89dd506456048be1cad58f6d383a9"
printf 'NGA_HOME=%s\n' "$dest"
