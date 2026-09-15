#!/usr/bin/env python3
from pathlib import Path


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected 1 match, found {count}")
    return text.replace(old, new, 1)


hpath = Path("B0Trackers.h")
h = hpath.read_text()
h = replace_once(
    h,
    "    bool m_failOnIncompleteSurfaceMap = false;\n",
    "    bool m_failOnIncompleteSurfaceMap = false;\n"
    "    bool m_enableTruthSeededChain = true;\n"
    "    bool m_enableStubSeededChain = true;\n"
    "    bool m_writeTrackStates = true;\n",
    "runtime config members",
)
hpath.write_text(h)

ccpath = Path("B0Trackers.cc")
cc = ccpath.read_text()
cc = replace_once(
    cc,
    "    app->SetDefaultParameter(\"B0Trackers:fail_on_incomplete_surface_map\",\n                             m_failOnIncompleteSurfaceMap,\n                             \"Throw in Init unless every B0 sensor maps to an ACTS surface\");\n",
    "    app->SetDefaultParameter(\"B0Trackers:fail_on_incomplete_surface_map\",\n                             m_failOnIncompleteSurfaceMap,\n                             \"Throw in Init unless every B0 sensor maps to an ACTS surface\");\n"
    "    app->SetDefaultParameter(\"B0Trackers:enable_truth_seeded_chain\", m_enableTruthSeededChain,\n"
    "                             \"Request and analyze the truth-seeded B0 CKF chain\");\n"
    "    app->SetDefaultParameter(\"B0Trackers:enable_stub_seeded_chain\", m_enableStubSeededChain,\n"
    "                             \"Request and analyze the stub-seeded B0 CKF chain\");\n"
    "    app->SetDefaultParameter(\"B0Trackers:write_track_states\", m_writeTrackStates,\n"
    "                             \"Request ACTS track containers and write per-state diagnostics\");\n",
    "runtime parameters",
)
cc = replace_once(
    cc,
    "    m_tree->Branch(\"eventNumber\", &m_eventNumber);\n",
    "    m_tree->Branch(\"eventNumber\", &m_eventNumber);\n"
    "    m_tree->Branch(\"config_enable_truth_seeded_chain\", &m_enableTruthSeededChain);\n"
    "    m_tree->Branch(\"config_enable_stub_seeded_chain\", &m_enableStubSeededChain);\n"
    "    m_tree->Branch(\"config_write_track_states\", &m_writeTrackStates);\n",
    "runtime metadata branches",
)

cc = replace_once(
    cc,
    "    const bool hasStubSeeds  = getOpt(event, \"B0TrackerSeeds\", stubSeeds);\n    const bool hasTruthSeeds = getOpt(event, \"B0TrackerTruthSeeds\", truthSeeds);\n",
    "    bool hasStubSeeds = false;\n"
    "    bool hasTruthSeeds = false;\n"
    "    if (m_enableStubSeededChain) hasStubSeeds = getOpt(event, \"B0TrackerSeeds\", stubSeeds);\n"
    "    if (m_enableTruthSeededChain) hasTruthSeeds = getOpt(event, \"B0TrackerTruthSeeds\", truthSeeds);\n",
    "conditional seed gets",
)

old_ts = """    const bool hasTsTrackParams =
        getOpt(event, "B0TrackerCKFTruthSeededTrackParameters", tsTracks);
    const bool hasTsTrajectories =
        getOpt(event, "B0TrackerCKFTruthSeededTrajectories", tsTrajectories);
    const bool hasTsTrajectoriesUnfiltered =
        getOpt(event, "B0TrackerCKFTruthSeededTrajectoriesUnfiltered", tsTrajectoriesUnfiltered);
    const bool hasTsTracks = getOpt(event, "B0TrackerCKFTruthSeededTracks", tsEdmTracks);
    const bool hasTsAssocs =
        getOpt(event, "B0TrackerCKFTruthSeededTrackAssociations", tsAssocs);
    const bool hasTsActsStates =
        getOpt(event, "B0TrackerCKFTruthSeededActsTrackStates", tsActsTrackStates);
    const bool hasTsActsTracks =
        getOpt(event, "B0TrackerCKFTruthSeededActsTracks", tsActsTracks);
    const bool hasTsTracksUnfiltered =
        getOpt(event, "B0TrackerCKFTruthSeededTracksUnfiltered", tsUnfiltered);
"""
new_ts = """    bool hasTsTrackParams = false;
    bool hasTsTrajectories = false;
    bool hasTsTrajectoriesUnfiltered = false;
    bool hasTsTracks = false;
    bool hasTsAssocs = false;
    bool hasTsActsStates = false;
    bool hasTsActsTracks = false;
    bool hasTsTracksUnfiltered = false;
    if (m_enableTruthSeededChain) {
        hasTsTrackParams = getOpt(event, "B0TrackerCKFTruthSeededTrackParameters", tsTracks);
        hasTsTrajectories = getOpt(event, "B0TrackerCKFTruthSeededTrajectories", tsTrajectories);
        hasTsTrajectoriesUnfiltered = getOpt(
            event, "B0TrackerCKFTruthSeededTrajectoriesUnfiltered", tsTrajectoriesUnfiltered);
        hasTsTracks = getOpt(event, "B0TrackerCKFTruthSeededTracks", tsEdmTracks);
        hasTsAssocs = getOpt(event, "B0TrackerCKFTruthSeededTrackAssociations", tsAssocs);
        hasTsTracksUnfiltered = getOpt(event, "B0TrackerCKFTruthSeededTracksUnfiltered", tsUnfiltered);
        if (m_writeTrackStates) {
            hasTsActsStates = getOpt(event, "B0TrackerCKFTruthSeededActsTrackStates", tsActsTrackStates);
            hasTsActsTracks = getOpt(event, "B0TrackerCKFTruthSeededActsTracks", tsActsTracks);
        }
    }
"""
cc = replace_once(cc, old_ts, new_ts, "conditional truth-chain gets")

old_ckf = """    const bool hasCkfTrackParams =
        getOpt(event, "B0TrackerCKFTrackParameters", ckfTracks);
    const bool hasCkfTrajectories = getOpt(event, "B0TrackerCKFTrajectories", ckfTrajectories);
    const bool hasCkfTrajectoriesUnfiltered =
        getOpt(event, "B0TrackerCKFTrajectoriesUnfiltered", ckfTrajectoriesUnfiltered);
    const bool hasCkfTracks  = getOpt(event, "B0TrackerCKFTracks", ckfEdmTracks);
    const bool hasCkfAssocs  = getOpt(event, "B0TrackerCKFTrackAssociations", ckfAssocs);
    const bool hasCkfActsStates =
        getOpt(event, "B0TrackerCKFActsTrackStates", ckfActsTrackStates);
    const bool hasCkfActsTracks = getOpt(event, "B0TrackerCKFActsTracks", ckfActsTracks);
    const bool hasCkfTracksUnfiltered =
        getOpt(event, "B0TrackerCKFTracksUnfiltered", ckfUnfiltered);
"""
new_ckf = """    bool hasCkfTrackParams = false;
    bool hasCkfTrajectories = false;
    bool hasCkfTrajectoriesUnfiltered = false;
    bool hasCkfTracks = false;
    bool hasCkfAssocs = false;
    bool hasCkfActsStates = false;
    bool hasCkfActsTracks = false;
    bool hasCkfTracksUnfiltered = false;
    if (m_enableStubSeededChain) {
        hasCkfTrackParams = getOpt(event, "B0TrackerCKFTrackParameters", ckfTracks);
        hasCkfTrajectories = getOpt(event, "B0TrackerCKFTrajectories", ckfTrajectories);
        hasCkfTrajectoriesUnfiltered = getOpt(
            event, "B0TrackerCKFTrajectoriesUnfiltered", ckfTrajectoriesUnfiltered);
        hasCkfTracks = getOpt(event, "B0TrackerCKFTracks", ckfEdmTracks);
        hasCkfAssocs = getOpt(event, "B0TrackerCKFTrackAssociations", ckfAssocs);
        hasCkfTracksUnfiltered = getOpt(event, "B0TrackerCKFTracksUnfiltered", ckfUnfiltered);
        if (m_writeTrackStates) {
            hasCkfActsStates = getOpt(event, "B0TrackerCKFActsTrackStates", ckfActsTrackStates);
            hasCkfActsTracks = getOpt(event, "B0TrackerCKFActsTracks", ckfActsTracks);
        }
    }
"""
cc = replace_once(cc, old_ckf, new_ckf, "conditional stub-chain gets")
ccpath.write_text(cc)
