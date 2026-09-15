from pathlib import Path

h = Path("B0Trackers.h")
text = h.read_text()
text = text.replace(
    "        int seedIndex = -1;\n        int identityValid = 0;",
    "        int seedIndex = -1;\n        std::uint32_t seedCollectionID = 0;\n        int identityValid = 0;",
    1,
)
text = text.replace(
    "        std::vector<int>    object_index, seed_index, identity_valid;\n        std::vector<std::uint32_t> object_collectionID;",
    "        std::vector<int>    object_index, seed_index, identity_valid;\n        std::vector<std::uint32_t> object_collectionID, seed_collectionID;",
    1,
)
text = text.replace(
    "        std::vector<int>    state_parent_seed_index, state_parent_track_index, state_parent_identity_valid;\n        std::vector<std::uint32_t> state_parent_track_collectionID;",
    "        std::vector<int>    state_parent_seed_index, state_parent_track_index, state_parent_identity_valid;\n        std::vector<std::uint32_t> state_parent_seed_collectionID, state_parent_track_collectionID;",
    1,
)
h.write_text(text)

p = Path("B0Trackers.cc")
text = p.read_text()
text = text.replace(
    "        std::map<int, std::vector<std::pair<std::uint32_t, int>>> trackObjectsBySeed;",
    "        std::map<std::pair<std::uint32_t, int>, std::vector<std::pair<std::uint32_t, int>>> trackObjectsBySeed;",
    1,
)
text = text.replace(
    "            if (seed.isAvailable()) {\n                trackObjectsBySeed[seed.id().index].push_back({objectId.collectionID, objectId.index});\n            }",
    "            if (seed.isAvailable()) {\n                const auto seedId = seed.id();\n                trackObjectsBySeed[{seedId.collectionID, seedId.index}].push_back(\n                    {objectId.collectionID, objectId.index});\n            }",
    1,
)
text = text.replace(
    "                         const auto& assocs,\n                         const auto& actsTracks,\n                         const auto& actsStates) {",
    "                         const auto& assocs,\n                         const auto& actsTracks,\n                         const auto& actsStates,\n                         const auto& chainSeeds) {",
    1,
)
text = text.replace(
    "            int seedIndex = -1;\n            int identityValid = 0;",
    "            int seedIndex = -1;\n            std::uint32_t seedCollectionID = 0;\n            int identityValid = 0;",
    1,
)
text = text.replace(
    "            const auto trajectorySeed = trajectory->getSeed();\n            if (trajectorySeed.isAvailable()) seedIndex = trajectorySeed.id().index;",
    "            const auto trajectorySeed = trajectory->getSeed();\n            if (trajectorySeed.isAvailable()) {\n                const auto seedId = trajectorySeed.id();\n                seedIndex = seedId.index;\n                seedCollectionID = seedId.collectionID;\n            }",
    1,
)
text = text.replace(
    "            out.seed_index.push_back(seedIndex);\n            out.identity_valid.push_back(identityValid);",
    "            out.seed_index.push_back(seedIndex);\n            out.seed_collectionID.push_back(seedCollectionID);\n            out.identity_valid.push_back(identityValid);",
    1,
)
text = text.replace(
    "                out.oracle.seedIndex = seedIndex; out.oracle.identityValid = identityValid;",
    "                out.oracle.seedIndex = seedIndex; out.oracle.seedCollectionID = seedCollectionID;\n                out.oracle.identityValid = identityValid;",
    1,
)
text = text.replace(
    "                    out.truthMatched.seedIndex = seedIndex;\n                    out.truthMatched.identityValid = identityValid;",
    "                    out.truthMatched.seedIndex = seedIndex;\n                    out.truthMatched.seedCollectionID = seedCollectionID;\n                    out.truthMatched.identityValid = identityValid;",
    1,
)
text = text.replace(
    "                    out.recoBest.seedIndex = seedIndex;\n                    out.recoBest.identityValid = identityValid;",
    "                    out.recoBest.seedIndex = seedIndex;\n                    out.recoBest.seedCollectionID = seedCollectionID;\n                    out.recoBest.identityValid = identityValid;",
    1,
)
text = text.replace(
    "            int parentSeedIndex = -1;\n            try { parentSeedIndex = static_cast<int>(seedNumber(track)); } catch (...) {}\n            int parentTrackIndex = -1;",
    "            int parentSeedIndex = -1;\n            std::uint32_t parentSeedCollectionID = 0;\n            try {\n                const auto seedPosition = static_cast<std::size_t>(seedNumber(track));\n                if (seedPosition < chainSeeds.size() && chainSeeds[seedPosition] != nullptr) {\n                    const auto seedId = chainSeeds[seedPosition]->id();\n                    parentSeedIndex = seedId.index;\n                    parentSeedCollectionID = seedId.collectionID;\n                }\n            } catch (...) {}\n            int parentTrackIndex = -1;",
    1,
)
text = text.replace(
    "            if (const auto found = trackObjectsBySeed.find(parentSeedIndex);\n                found != trackObjectsBySeed.end() && found->second.size() == 1) {",
    "            if (const auto found = trackObjectsBySeed.find({parentSeedCollectionID, parentSeedIndex});\n                parentSeedIndex >= 0 && found != trackObjectsBySeed.end() &&\n                found->second.size() == 1) {",
    1,
)
text = text.replace(
    "                out.state_parent_seed_index.push_back(parentSeedIndex);\n                out.state_parent_track_index.push_back(parentTrackIndex);",
    "                out.state_parent_seed_index.push_back(parentSeedIndex);\n                out.state_parent_seed_collectionID.push_back(parentSeedCollectionID);\n                out.state_parent_track_index.push_back(parentTrackIndex);",
    1,
)
text = text.replace(
    "    fillChain(m_ts, tsTrajectories, tsTracks, tsEdmTracks, tsAssocs, tsActsTracks, tsActsTrackStates);\n    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks, ckfActsTrackStates);",
    "    fillChain(m_ts, tsTrajectories, tsTracks, tsEdmTracks, tsAssocs, tsActsTracks,\n              tsActsTrackStates, truthSeeds);\n    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks,\n              ckfActsTrackStates, stubSeeds);",
    1,
)
text = text.replace(
    "    objectIndex = -1; objectCollectionID = 0; seedIndex = -1; identityValid = 0;",
    "    objectIndex = -1; objectCollectionID = 0; seedIndex = -1; seedCollectionID = 0; identityValid = 0;",
    1,
)
text = text.replace(
    "    object_index.clear(); object_collectionID.clear(); seed_index.clear(); identity_valid.clear();",
    "    object_index.clear(); object_collectionID.clear(); seed_index.clear(); seed_collectionID.clear();\n    identity_valid.clear();",
    1,
)
text = text.replace(
    "    state_parent_seed_index.clear(); state_parent_track_index.clear();",
    "    state_parent_seed_index.clear(); state_parent_seed_collectionID.clear();\n    state_parent_track_index.clear();",
    1,
)
text = text.replace(
    "    br(prefix + \"seed_index\", &b.seedIndex);\n    br(prefix + \"identity_valid\", &b.identityValid);",
    "    br(prefix + \"seed_index\", &b.seedIndex);\n    br(prefix + \"seed_collectionID\", &b.seedCollectionID);\n    br(prefix + \"identity_valid\", &b.identityValid);",
    1,
)
text = text.replace(
    "    br(trkPrefix + \"seed_index\", &c.seed_index);\n    br(trkPrefix + \"identity_valid\", &c.identity_valid);",
    "    br(trkPrefix + \"seed_index\", &c.seed_index);\n    br(trkPrefix + \"seed_collectionID\", &c.seed_collectionID);\n    br(trkPrefix + \"identity_valid\", &c.identity_valid);",
    1,
)
text = text.replace(
    "    br(trkPrefix + \"state_parent_seed_index\", &c.state_parent_seed_index);\n    br(trkPrefix + \"state_parent_track_index\", &c.state_parent_track_index);",
    "    br(trkPrefix + \"state_parent_seed_index\", &c.state_parent_seed_index);\n    br(trkPrefix + \"state_parent_seed_collectionID\", &c.state_parent_seed_collectionID);\n    br(trkPrefix + \"state_parent_track_index\", &c.state_parent_track_index);",
    1,
)
p.write_text(text)
