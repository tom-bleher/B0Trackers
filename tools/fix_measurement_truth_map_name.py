from pathlib import Path

path = Path("B0Trackers.cc")
text = path.read_text()
old = "const auto truthIt = cellTruthEDep.find(cid);\n                if (truthIt == cellTruthEDep.end()) continue;"
new = "const auto truthIt = cellParticleEDepByCell.find(cid);\n                if (truthIt == cellParticleEDepByCell.end()) continue;"
if old not in text:
    raise SystemExit("measurement truth map reference not found")
path.write_text(text.replace(old, new, 1))
