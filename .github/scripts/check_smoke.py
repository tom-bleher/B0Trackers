"""Assert that a B0Trackers smoke-run tree is complete and self-consistent."""

import re
import sys
from pathlib import Path

import awkward as ak
import uproot

# The schema version the plugin was built with lives in the C++ helper header.
HEADER = Path(__file__).resolve().parents[2] / "B0TrackersHelpers.h"
SCHEMA_VERSION = int(re.search(r"kSchemaVersion\s*=\s*(\d+)", HEADER.read_text()).group(1))

tree = uproot.open(sys.argv[1])["B0Trackers/hits"]
failures = []

n = tree.num_entries
if n != 5:
    failures.append(f"expected 5 entries, got {n}")

schema = set(tree["schema_version"].array(library="np").tolist())
if schema != {SCHEMA_VERSION}:
    failures.append(f"schema_version {schema}, expected {{{SCHEMA_VERSION}}}")

# Every collection the stock EICrecon registers must be reported present: a
# false here is how the JANA2 GetFactory/GetCollection change hid itself.
for name in [k for k in tree.keys() if k.startswith("has_")]:
    if not ak.all(tree[name].array()):
        failures.append(f"{name} is false")

if int(ak.sum(tree["n_simhits"].array())) == 0:
    failures.append("no B0 SimTrackerHits: the gun missed the B0 tracker")

station = ak.flatten(tree["station"].array(), axis=None)
max_station = int(ak.max(station)) if len(station) else 0
if not 1 <= max_station <= 4:
    failures.append(f"station range ends at {max_station}, expected 1-4")

for failure in failures:
    print("SMOKE FAIL:", failure)
print(f"smoke: {n} events, schema {sorted(schema)}, {int(ak.sum(tree['n_simhits'].array()))} B0 SimHits,"
      f" stations up to {max_station}")
sys.exit(1 if failures else 0)
