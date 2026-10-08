# Sourced by the pipeline shell scripts: workspace paths.
#   KANSASII_ROOT  workspace root (override e.g. for a mapped drive)
#   ARM_DIR        where refs/ and work/ live (default: $KANSASII_ROOT/output/arm)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}
ARM=${ARM_DIR:-$ROOT/output/arm}
