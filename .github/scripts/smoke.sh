#!/bin/bash
# Run the plugin on a few B0 protons in the eic-shell stock geometry and
# EICrecon, and check the tree: the build alone does not catch a JANA2 or
# geometry-naming change that only shows up at runtime.
#   .github/scripts/smoke.sh <plugin-install-dir> <work-dir>
set -euo pipefail
PLUGINS=$(realpath "$1")
WORK=$(realpath -m "$2")
HERE=$(dirname "$(realpath "$0")")
mkdir -p "$WORK"
cd "$WORK"

: "${DETECTOR_PATH:?DETECTOR_PATH is not set; run inside eic-shell}"
# A non-interactive eic-shell does not load the geometry plugins; source the
# thisepic.sh that belongs to DETECTOR_PATH (<prefix>/share/epic).
THISEPIC=${THISEPIC:-$(dirname "$(dirname "$DETECTOR_PATH")")/bin/thisepic.sh}
if [ -f "$THISEPIC" ]; then
  # shellcheck disable=SC1090
  source "$THISEPIC"
fi
# Pinned rather than taken from DETECTOR_CONFIG, which eic-shell does not
# always set to a geometry.
CONFIG=${B0_SMOKE_CONFIG:-epic_craterlake}
test -f "$DETECTOR_PATH/$CONFIG.xml" || { echo "no $DETECTOR_PATH/$CONFIG.xml"; exit 1; }
export DETECTOR_CONFIG=$CONFIG
echo "geometry: $DETECTOR_PATH/$CONFIG.xml"

# The ion beam runs at -25 mrad in x, so aim the gun just off it (phi = pi).
npsim --compactFile "$DETECTOR_PATH/$CONFIG.xml" \
  --enableGun --gun.particle proton \
  --gun.momentumMin "80*GeV" --gun.momentumMax "100*GeV" \
  --gun.thetaMin "0.030*rad" --gun.thetaMax "0.040*rad" \
  --gun.phiMin "175*deg" --gun.phiMax "185*deg" \
  --gun.distribution uniform \
  --numberOfEvents 5 --outputFile sim.edm4hep.root > npsim.log 2>&1 \
  || { tail -30 npsim.log; exit 1; }

JANA_PLUGIN_PATH="$PLUGINS${JANA_PLUGIN_PATH:+:$JANA_PLUGIN_PATH}" \
eicrecon -Pplugins=B0Trackers -Pnthreads=2 \
  -Phistsfile=b0trackers.root -Ppodio:output_file=rec.edm4eic.root \
  sim.edm4hep.root > eicrecon.log 2>&1 \
  || { grep -B2 -A8 'JException\|fatal' eicrecon.log | tail -40; exit 1; }
grep '\[B0Trackers\]' eicrecon.log | head -5

python3 "$HERE/check_smoke.py" b0trackers.root
