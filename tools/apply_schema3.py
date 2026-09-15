#!/usr/bin/env python3
from pathlib import Path


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {count}")
    return text.replace(old, new, 1)


hpath = Path("B0Trackers.h")
h = hpath.read_text()
h = replace_once(h, "#include <cstdint>\n", "#include <array>\n#include <cstdint>\n", "header array include")
h = replace_once(
    h,
    "        int assocMcIndex = -1;\n        int pdg = 0;",
    "        int assocMcIndex = -1;\n        std::uint32_t assocMcCollectionID = 0;\n        int objectIndex = -1;\n        std::uint32_t objectCollectionID = 0;\n        int seedIndex = -1;\n        int identityValid = 0;\n        int pdg = 0;",
    "BestSel identity",
)
h = replace_once(
    h,
    "        std::vector<int>    index, charge, type, pdg;\n        std::vector<int>    nStates, nMeasurements, nOutliers, nHoles, nSharedHits;",
    "        std::vector<int>    index, charge, type, pdg;\n        std::vector<int>    object_index, seed_index, identity_valid;\n        std::vector<std::uint32_t> object_collectionID;\n        std::array<std::vector<double>, 21> cov_upper;\n        std::vector<int>    nStates, nMeasurements, nOutliers, nHoles, nSharedHits;",
    "track identity/covariance",
)
h = replace_once(
    h,
    "        std::vector<int>    state_track_index, state_index, state_acts_index, state_type, state_pdg;\n        std::vector<int>    state_mapping_method;",
    "        std::vector<int>    state_track_index, state_index, state_acts_index, state_type, state_pdg;\n        std::vector<int>    state_parent_seed_index, state_parent_track_index, state_parent_identity_valid;\n        std::vector<std::uint32_t> state_parent_track_collectionID;\n        std::vector<int>    state_estimate_kind;\n        std::vector<int>    state_mapping_method;",
    "state identity",
)
h = replace_once(
    h,
    "        std::vector<double> state_meas_loc0, state_meas_loc1, state_resid_loc0, state_resid_loc1;\n        std::vector<double> x_on_plane, y_on_plane, z_on_plane;",
    "        std::vector<double> state_meas_loc0, state_meas_loc1, state_resid_loc0, state_resid_loc1;\n        std::vector<int>    state_meas_dim, state_proj_index0, state_proj_index1;\n        std::vector<double> state_meas_cov00, state_meas_cov01, state_meas_cov11;\n        std::vector<double> state_pred_loc0, state_pred_loc1, state_pred_theta, state_pred_phi;\n        std::vector<double> state_pred_qOverP, state_pred_time;\n        std::array<std::vector<double>, 21> state_pred_cov_upper;\n        std::vector<double> state_innov0, state_innov1, state_innov_pull0, state_innov_pull1;\n        std::vector<double> state_innov_chi2;\n        std::vector<double> x_on_plane, y_on_plane, z_on_plane;",
    "innovation vectors",
)
h = replace_once(
    h,
    "        std::vector<double> state_theta, state_phi, state_qOverP, state_time;\n\n        BestSel oracle;",
    "        std::vector<double> state_theta, state_phi, state_qOverP, state_time;\n\n        int nMapExact = 0;\n        int nMapFallback = 0;\n        int nMapFailedPhysics = 0;\n        int nMapUnmappedNonPhysics = 0;\n\n        BestSel oracle;",
    "per-chain mapping counters",
)
h = replace_once(
    h,
    "    std::vector<int>    vm_seed_n_unfiltered_tracks, vm_seed_n_filtered_tracks;\n    std::vector<double> vm_truth_seed_quality",
    "    std::vector<int>    vm_seed_n_unfiltered_tracks, vm_seed_n_filtered_tracks;\n    std::vector<int>    vm_seed_assoc_mcIndex;\n    std::vector<std::uint32_t> vm_seed_assoc_mcCollectionID;\n    std::vector<double> vm_seed_assoc_weight;\n    std::vector<double> vm_truth_seed_quality",
    "seed truth branches",
)
h = replace_once(
    h,
    "    int m_nStationsPrimary = 0; // selected primary only\n\n    int m_nSimHits = 0;",
    "    int m_nStationsPrimary = 0; // selected primary only\n    int m_selPrimaryHasSeed = 0;\n    int m_selPrimaryHasUnfilteredTrack = 0;\n    int m_selPrimaryHasFilteredTrack = 0;\n    int m_selPrimaryHasTruthMatchedTrack = 0;\n\n    int m_nSimHits = 0;",
    "selected-primary stages",
)
h = replace_once(
    h,
    "    bool m_hasCkfTracksUnfiltered = false;\n};",
    "    bool m_hasCkfTracksUnfiltered = false;\n    bool m_hasCkfAssocsUnfiltered = false;\n};",
    "unfiltered assoc availability",
)
hpath.write_text(h)

ccpath = Path("B0Trackers.cc")
cc = ccpath.read_text()
cc = replace_once(cc, "#include <algorithm>\n", "#include <algorithm>\n#include <array>\n", "cc array include")
cc = replace_once(cc, "#include <utility>\n", "#include <utility>\n\n#include <Eigen/Cholesky>\n", "eigen include")
cc = replace_once(
    cc,
    "#include <Acts/Definitions/TrackParametrization.hpp>\n",
    "#include <Acts/Definitions/TrackParametrization.hpp>\n#include <Acts/Definitions/Units.hpp>\n",
    "acts units include",
)
cc = replace_once(
    cc,
    "#include <Acts/EventData/TrackContainer.hpp>\n",
    "#include <Acts/EventData/ProxyAccessor.hpp>\n#include <Acts/EventData/TrackContainer.hpp>\n",
    "proxy accessor include",
)
cc = cc.replace("// B0Trackers/hits branch map (schema 2):", "// B0Trackers/hits branch map (schema 3):", 1)

cc = replace_once(
    cc,
    "    m_tree->Branch(\"n_stations_primary\", &m_nStationsPrimary);\n",
    "    m_tree->Branch(\"n_stations_primary\", &m_nStationsPrimary);\n"
    "    m_tree->Branch(\"sel_primary_has_seed\", &m_selPrimaryHasSeed);\n"
    "    m_tree->Branch(\"sel_primary_has_unfiltered_track\", &m_selPrimaryHasUnfilteredTrack);\n"
    "    m_tree->Branch(\"sel_primary_has_filtered_track\", &m_selPrimaryHasFilteredTrack);\n"
    "    m_tree->Branch(\"sel_primary_has_truth_matched_track\", &m_selPrimaryHasTruthMatchedTrack);\n",
    "selected stage branches",
)
cc = replace_once(
    cc,
    "    m_tree->Branch(\"seed_n_filtered_tracks\", &vm_seed_n_filtered_tracks);\n",
    "    m_tree->Branch(\"seed_n_filtered_tracks\", &vm_seed_n_filtered_tracks);\n"
    "    m_tree->Branch(\"seed_assoc_mcIndex\", &vm_seed_assoc_mcIndex);\n"
    "    m_tree->Branch(\"seed_assoc_mcCollectionID\", &vm_seed_assoc_mcCollectionID);\n"
    "    m_tree->Branch(\"seed_assoc_weight\", &vm_seed_assoc_weight);\n",
    "seed assoc branches",
)
cc = replace_once(
    cc,
    "    m_tree->Branch(\"has_ckf_tracks_unfiltered\", &m_hasCkfTracksUnfiltered);\n",
    "    m_tree->Branch(\"has_ckf_tracks_unfiltered\", &m_hasCkfTracksUnfiltered);\n"
    "    m_tree->Branch(\"has_ckf_assocs_unfiltered\", &m_hasCkfAssocsUnfiltered);\n",
    "unfiltered assoc branch",
)

cc = replace_once(
    cc,
    "    std::vector<const edm4eic::Track*> ckfUnfiltered;\n",
    "    std::vector<const edm4eic::Track*> ckfUnfiltered;\n"
    "    std::vector<const edm4eic::MCRecoTrackParticleAssociation*> ckfUnfilteredAssocs;\n",
    "unfiltered assoc vector",
)
cc = replace_once(
    cc,
    "    const bool hasCkfTracksUnfiltered =\n        getOpt(event, \"B0TrackerCKFTracksUnfiltered\", ckfUnfiltered);\n",
    "    const bool hasCkfTracksUnfiltered =\n        getOpt(event, \"B0TrackerCKFTracksUnfiltered\", ckfUnfiltered);\n"
    "    const bool hasCkfAssocsUnfiltered =\n        getOpt(event, \"B0TrackerCKFTrackUnfilteredAssociations\", ckfUnfilteredAssocs);\n",
    "unfiltered assoc get",
)
cc = replace_once(
    cc,
    "    m_hasCkfTracksUnfiltered     = hasCkfTracksUnfiltered;\n",
    "    m_hasCkfTracksUnfiltered     = hasCkfTracksUnfiltered;\n"
    "    m_hasCkfAssocsUnfiltered     = hasCkfAssocsUnfiltered;\n",
    "unfiltered assoc assignment",
)
cc = replace_once(
    cc,
    "    vm_seed_n_unfiltered_tracks.clear(); vm_seed_n_filtered_tracks.clear();\n",
    "    vm_seed_n_unfiltered_tracks.clear(); vm_seed_n_filtered_tracks.clear();\n"
    "    vm_seed_assoc_mcIndex.clear(); vm_seed_assoc_mcCollectionID.clear(); vm_seed_assoc_weight.clear();\n",
    "seed assoc clear",
)
cc = replace_once(
    cc,
    "    m_nStationsPrimary = 0;\n",
    "    m_nStationsPrimary = 0;\n"
    "    m_selPrimaryHasSeed = 0;\n"
    "    m_selPrimaryHasUnfilteredTrack = 0;\n"
    "    m_selPrimaryHasFilteredTrack = 0;\n"
    "    m_selPrimaryHasTruthMatchedTrack = 0;\n",
    "stage reset",
)

cc = replace_once(
    cc,
    "    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;\n    {\n        // (cellID, MC ObjectID) -> summed deposit.\n        std::map<std::pair<std::uint64_t, std::pair<std::uint32_t, int>>, double> cellParticleEDep;",
    "    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;\n"
    "    std::unordered_map<std::uint64_t, std::map<std::pair<std::uint32_t, int>, double>> cellParticleEDepByCell;\n"
    "    {\n        // (cellID, MC ObjectID) -> summed deposit.\n"
    "        std::map<std::pair<std::uint64_t, std::pair<std::uint32_t, int>>, double> cellParticleEDep;",
    "cell contribution map declaration",
)
cc = replace_once(
    cc,
    "            cellParticleEDep[{cid, {id.collectionID, id.index}}] += sim.getEDep();\n",
    "            cellParticleEDep[{cid, {id.collectionID, id.index}}] += sim.getEDep();\n"
    "            cellParticleEDepByCell[cid][{id.collectionID, id.index}] += sim.getEDep();\n",
    "cell contribution fill",
)

seed_anchor = """    fillSeeds(truthSeeds, tsTrajectories, tsTrajectoriesUnfiltered, hasTsTrajectoriesUnfiltered,
              vm_truth_seed_quality, vm_truth_seed_p, vm_truth_seed_qOverP,
              vm_truth_seed_theta, vm_truth_seed_phi, vm_truth_seed_loc0, vm_truth_seed_loc1,
              vm_truth_seed_sigma_qOverP, vm_truth_seed_sigma_theta, vm_truth_seed_sigma_phi,
              vm_truth_seed_nHits, vm_truth_seed_charge, vm_truth_seed_momentum_resolved,
              vm_truth_seed_became_track,
              vm_truth_seed_made_unfiltered_track, vm_truth_seed_survived_ambiguity,
              vm_truth_seed_n_unfiltered_tracks, vm_truth_seed_n_filtered_tracks);

"""
seed_insert = seed_anchor + """    // Attribute each stub seed to MC truth from its constituent TrackerHits.
    // Each hit/cell contributes one unit distributed among the cell's MC energy fractions,
    // matching the per-measurement convention used by ActsToTracks.
    for (const auto* seed : stubSeeds) {
        std::map<std::pair<std::uint32_t, int>, double> weights;
        if (seed != nullptr) {
            for (std::size_t ih = 0; ih < seed->hits_size(); ++ih) {
                const auto hit = seed->getHits(ih);
                if (!hit.isAvailable()) continue;
                const auto found = cellParticleEDepByCell.find(hit.getCellID());
                if (found == cellParticleEDepByCell.end()) continue;
                double total = 0.0;
                for (const auto& [id, edep] : found->second) total += edep;
                if (!(total > 0.0)) continue;
                for (const auto& [id, edep] : found->second) weights[id] += edep / total;
            }
        }
        int bestIndex = -1;
        std::uint32_t bestCollection = 0;
        double bestWeight = -1.0;
        double totalWeight = 0.0;
        for (const auto& [id, weight] : weights) {
            totalWeight += weight;
            if (weight > bestWeight) {
                bestWeight = weight;
                bestCollection = id.first;
                bestIndex = id.second;
            }
        }
        const double normalized = totalWeight > 0.0 ? bestWeight / totalWeight : nan;
        vm_seed_assoc_mcIndex.push_back(bestIndex);
        vm_seed_assoc_mcCollectionID.push_back(bestCollection);
        vm_seed_assoc_weight.push_back(normalized);
        if (bestIndex == m_selPrimaryMcIndex && bestCollection == m_selPrimaryMcCollectionID &&
            m_selPrimaryMcIndex >= 0) {
            m_selPrimaryHasSeed = 1;
        }
    }

    const auto collectionHasSelectedPrimary = [this](const auto& assocs) {
        std::map<std::pair<std::uint32_t, int>, std::tuple<double, std::uint32_t, int>> best;
        for (const auto* assoc : assocs) {
            if (assoc == nullptr) continue;
            const auto rec = assoc->getRec();
            const auto sim = assoc->getSim();
            if (!rec.isAvailable() || !sim.isAvailable()) continue;
            const auto rid = rec.id();
            const auto sid = sim.id();
            const auto key = std::make_pair(rid.collectionID, rid.index);
            const double weight = assoc->getWeight();
            const auto it = best.find(key);
            if (it == best.end() || weight > std::get<0>(it->second)) {
                best[key] = {weight, sid.collectionID, sid.index};
            }
        }
        for (const auto& [key, value] : best) {
            if (std::get<1>(value) == m_selPrimaryMcCollectionID &&
                std::get<2>(value) == m_selPrimaryMcIndex && m_selPrimaryMcIndex >= 0) return true;
        }
        return false;
    };
    m_selPrimaryHasUnfilteredTrack =
        hasCkfAssocsUnfiltered && collectionHasSelectedPrimary(ckfUnfilteredAssocs) ? 1 : 0;
    m_selPrimaryHasFilteredTrack =
        hasCkfAssocs && collectionHasSelectedPrimary(ckfAssocs) ? 1 : 0;

"""
cc = replace_once(cc, seed_anchor, seed_insert, "seed truth/stage insertion")

# Stable track identity: associations are keyed by complete PODIO ObjectID and
# EDM tracks are resolved through their Trajectory relation instead of collection order.
cc = replace_once(
    cc,
    "        std::unordered_map<int, AssocHit> bestAssoc;\n",
    "        std::map<std::pair<std::uint32_t, int>, AssocHit> bestAssoc;\n",
    "bestAssoc key",
)
cc = replace_once(
    cc,
    "            const int recIndex = rec.id().index;\n            const double weight = assoc->getWeight();\n            auto it = bestAssoc.find(recIndex);\n            if (it == bestAssoc.end() || weight > it->second.weight) {\n                bestAssoc[recIndex] = {sim.id().index, sim.id().collectionID, weight};\n            }\n",
    "            const auto recId = rec.id();\n            const double weight = assoc->getWeight();\n            const auto recKey = std::make_pair(recId.collectionID, recId.index);\n            auto it = bestAssoc.find(recKey);\n            if (it == bestAssoc.end() || weight > it->second.weight) {\n                bestAssoc[recKey] = {sim.id().index, sim.id().collectionID, weight};\n            }\n",
    "bestAssoc build",
)
identity_anchor = """        const double primaryP = m_selPrimaryP;
"""
identity_insert = """        std::map<std::pair<std::uint32_t, int>, const edm4eic::Track*> edmTrackByTrajectory;
        std::map<int, std::vector<std::pair<std::uint32_t, int>>> trackObjectsBySeed;
        std::map<std::pair<std::uint32_t, int>, int> pdgByTrackObject;
        for (const auto* edmTrack : edmTracks) {
            if (edmTrack == nullptr) continue;
            const auto objectId = edmTrack->id();
            pdgByTrackObject[{objectId.collectionID, objectId.index}] = edmTrack->getPdg();
            const auto trajectory = edmTrack->getTrajectory();
            if (!trajectory.isAvailable()) continue;
            const auto trajectoryId = trajectory.id();
            edmTrackByTrajectory[{trajectoryId.collectionID, trajectoryId.index}] = edmTrack;
            const auto seed = trajectory.getSeed();
            if (seed.isAvailable()) {
                trackObjectsBySeed[seed.id().index].push_back({objectId.collectionID, objectId.index});
            }
        }

""" + identity_anchor
cc = replace_once(cc, identity_anchor, identity_insert, "identity maps")

cc = replace_once(
    cc,
    "            const int currentTrackIndex = static_cast<int>(trajIndex);\n            const int nStates = static_cast<int>(trajectory->getNStates());",
    "            const int currentTrackIndex = static_cast<int>(trajIndex); // legacy positional index\n"
    "            int objectIndex = -1;\n"
    "            std::uint32_t objectCollectionID = 0;\n"
    "            int seedIndex = -1;\n"
    "            int identityValid = 0;\n"
    "            const auto trajectoryId = trajectory->id();\n"
    "            const auto stableTrackIt = edmTrackByTrajectory.find({trajectoryId.collectionID, trajectoryId.index});\n"
    "            const edm4eic::Track* stableEdmTrack =\n"
    "                stableTrackIt == edmTrackByTrajectory.end() ? nullptr : stableTrackIt->second;\n"
    "            if (stableEdmTrack != nullptr) {\n"
    "                const auto objectId = stableEdmTrack->id();\n"
    "                objectIndex = objectId.index;\n"
    "                objectCollectionID = objectId.collectionID;\n"
    "                identityValid = 1;\n"
    "            }\n"
    "            const auto trajectorySeed = trajectory->getSeed();\n"
    "            if (trajectorySeed.isAvailable()) seedIndex = trajectorySeed.id().index;\n"
    "            const int nStates = static_cast<int>(trajectory->getNStates());",
    "track identity resolve",
)
cc = replace_once(
    cc,
    "            if (trajIndex < edmTracks.size() && edmTracks[trajIndex] != nullptr) {\n                chi2 = edmTracks[trajIndex]->getChi2();\n                ndf = static_cast<int>(edmTracks[trajIndex]->getNdf());\n                pdg = edmTracks[trajIndex]->getPdg();\n            }",
    "            if (stableEdmTrack != nullptr) {\n                chi2 = stableEdmTrack->getChi2();\n                ndf = static_cast<int>(stableEdmTrack->getNdf());\n                pdg = stableEdmTrack->getPdg();\n            } else if (trajIndex < edmTracks.size() && edmTracks[trajIndex] != nullptr) {\n                chi2 = edmTracks[trajIndex]->getChi2();\n                ndf = static_cast<int>(edmTracks[trajIndex]->getNdf());\n                pdg = edmTracks[trajIndex]->getPdg();\n                const auto objectId = edmTracks[trajIndex]->id();\n                objectIndex = objectId.index;\n                objectCollectionID = objectId.collectionID;\n            }",
    "edm stable use",
)
cc = replace_once(
    cc,
    "            if (const auto it = bestAssoc.find(currentTrackIndex); it != bestAssoc.end()) {\n                assocMc = it->second.mcIndex;\n                assocCol = it->second.mcCollectionID;\n                assocW = it->second.weight;\n            }",
    "            if (objectIndex >= 0) {\n                if (const auto it = bestAssoc.find({objectCollectionID, objectIndex}); it != bestAssoc.end()) {\n                    assocMc = it->second.mcIndex;\n                    assocCol = it->second.mcCollectionID;\n                    assocW = it->second.weight;\n                }\n            }",
    "association stable lookup",
)
cc = replace_once(
    cc,
    "            out.index.push_back(currentTrackIndex);\n            out.type.push_back(tp.getType());",
    "            out.index.push_back(currentTrackIndex);\n"
    "            out.object_index.push_back(objectIndex);\n"
    "            out.object_collectionID.push_back(objectCollectionID);\n"
    "            out.seed_index.push_back(seedIndex);\n"
    "            out.identity_valid.push_back(identityValid);\n"
    "            {\n"
    "                int k = 0;\n"
    "                for (unsigned i = 0; i < 6; ++i) {\n"
    "                    for (unsigned j = i; j < 6; ++j) out.cov_upper[k++].push_back(cov(i, j));\n"
    "                }\n"
    "            }\n"
    "            out.type.push_back(tp.getType());",
    "track identity push",
)
# Attach stable identity to each BestSel after the existing assignment.
cc = cc.replace(
    "                           assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);\n            }",
    "                           assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);\n"
    "                out.oracle.objectIndex = objectIndex; out.oracle.objectCollectionID = objectCollectionID;\n"
    "                out.oracle.seedIndex = seedIndex; out.oracle.identityValid = identityValid;\n"
    "                out.oracle.assocMcCollectionID = assocCol;\n"
    "            }",
    1,
)
truth_old = """                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                }
            }
            if (nMeasurements > 0) {
"""
truth_new = """                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                    out.truthMatched.objectIndex = objectIndex;
                    out.truthMatched.objectCollectionID = objectCollectionID;
                    out.truthMatched.seedIndex = seedIndex;
                    out.truthMatched.identityValid = identityValid;
                    out.truthMatched.assocMcCollectionID = assocCol;
                }
            }
            if (nMeasurements > 0) {
"""
cc = replace_once(cc, truth_old, truth_new, "truth best identity")
reco_old = """                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                }
            }
        }
"""
reco_new = """                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                    out.recoBest.objectIndex = objectIndex;
                    out.recoBest.objectCollectionID = objectCollectionID;
                    out.recoBest.seedIndex = seedIndex;
                    out.recoBest.identityValid = identityValid;
                    out.recoBest.assocMcCollectionID = assocCol;
                }
            }
        }
"""
cc = replace_once(cc, reco_old, reco_new, "reco best identity")

cc = replace_once(
    cc,
    "        const auto nActsTracks = static_cast<int>(trackContainer.size());\n        for (int actsTrackIndex = 0; actsTrackIndex < nActsTracks; ++actsTrackIndex) {",
    "        const auto nActsTracks = static_cast<int>(trackContainer.size());\n"
    "        Acts::ConstProxyAccessor<unsigned int> seedNumber(\"seed\");\n"
    "        for (int actsTrackIndex = 0; actsTrackIndex < nActsTracks; ++actsTrackIndex) {",
    "state seed accessor",
)
cc = replace_once(
    cc,
    "            const auto track = trackContainer.getTrack(actsTrackIndex);\n            const std::size_t stateBegin = out.state_index.size();",
    "            const auto track = trackContainer.getTrack(actsTrackIndex);\n"
    "            int parentSeedIndex = -1;\n"
    "            try { parentSeedIndex = static_cast<int>(seedNumber(track)); } catch (...) {}\n"
    "            int parentTrackIndex = -1;\n"
    "            std::uint32_t parentTrackCollectionID = 0;\n"
    "            int parentIdentityValid = 0;\n"
    "            if (const auto found = trackObjectsBySeed.find(parentSeedIndex);\n"
    "                found != trackObjectsBySeed.end() && found->second.size() == 1) {\n"
    "                parentTrackCollectionID = found->second.front().first;\n"
    "                parentTrackIndex = found->second.front().second;\n"
    "                parentIdentityValid = 1;\n"
    "            }\n"
    "            const std::size_t stateBegin = out.state_index.size();",
    "state parent identity",
)

estimate_old = """                const auto selectParams = [&state](auto&& fill) {
                    if (state.hasSmoothed()) {
                        fill(state.smoothed());
                    } else if (state.hasFiltered()) {
                        fill(state.filtered());
                    } else if (state.hasPredicted()) {
                        fill(state.predicted());
                    }
                };

                double loc0 = nan;
"""
estimate_new = """                const int estimateKind = b0trk::preferredEstimateKind(
                    state.hasPredicted(), state.hasFiltered(), state.hasSmoothed());
                const auto selectParams = [&state](auto&& fill) {
                    if (state.hasSmoothed()) {
                        fill(state.smoothed());
                    } else if (state.hasFiltered()) {
                        fill(state.filtered());
                    } else if (state.hasPredicted()) {
                        fill(state.predicted());
                    }
                };

                double loc0 = nan;
"""
cc = replace_once(cc, estimate_old, estimate_new, "estimate kind")
cc = replace_once(
    cc,
    "                    stateTime = params[Acts::eBoundTime];\n                });\n",
    "                    stateTime = b0trk::nativeTimeToNs(params[Acts::eBoundTime], Acts::UnitConstants::ns);\n                });\n",
    "state time ns",
)

innov_anchor = """                const auto flags   = state.typeFlags();
"""
innov_code = """                double predLoc0 = nan, predLoc1 = nan, predTheta = nan, predPhi = nan;
                double predQOverP = nan, predTime = nan;
                std::array<double, 21> predCovUpper{};
                predCovUpper.fill(nan);
                int measDim = 0, proj0 = -1, proj1 = -1;
                double measCov00 = nan, measCov01 = nan, measCov11 = nan;
                double innov0 = nan, innov1 = nan, innovPull0 = nan, innovPull1 = nan, innovChi2 = nan;
                if (state.hasPredicted()) {
                    const auto pred = state.predicted();
                    const auto predCov = state.predictedCovariance();
                    predLoc0 = pred[Acts::eBoundLoc0]; predLoc1 = pred[Acts::eBoundLoc1];
                    predPhi = pred[Acts::eBoundPhi]; predTheta = pred[Acts::eBoundTheta];
                    predQOverP = pred[Acts::eBoundQOverP];
                    predTime = b0trk::nativeTimeToNs(pred[Acts::eBoundTime], Acts::UnitConstants::ns);
                    const std::array<double, 6> scale{
                        1.0 / Acts::UnitConstants::mm, 1.0 / Acts::UnitConstants::mm,
                        1.0 / Acts::UnitConstants::rad, 1.0 / Acts::UnitConstants::rad,
                        Acts::UnitConstants::GeV, 1.0 / Acts::UnitConstants::ns};
                    int k = 0;
                    for (unsigned i = 0; i < 6; ++i) {
                        for (unsigned j = i; j < 6; ++j) predCovUpper[k++] = predCov(i, j) * scale[i] * scale[j];
                    }
                    if (state.hasCalibrated()) {
                        try {
                            const auto meas = state.effectiveCalibrated();
                            const auto measCov = state.effectiveCalibratedCovariance();
                            const auto projector = state.projectorSubspaceIndices();
                            measDim = static_cast<int>(state.calibratedSize());
                            if (measDim > 0) {
                                proj0 = static_cast<int>(projector[0]);
                                measCov00 = measCov(0, 0);
                                innov0 = meas[0] - pred[projector[0]];
                            }
                            if (measDim > 1) {
                                proj1 = static_cast<int>(projector[1]);
                                measCov01 = measCov(0, 1);
                                measCov11 = measCov(1, 1);
                                innov1 = meas[1] - pred[projector[1]];
                            }
                            if (measDim == 1) {
                                const double s00 = measCov(0, 0) + predCov(projector[0], projector[0]);
                                if (s00 > 0.0 && std::isfinite(s00)) {
                                    innovPull0 = innov0 / std::sqrt(s00);
                                    innovChi2 = innov0 * innov0 / s00;
                                }
                            } else if (measDim == 2) {
                                Eigen::Matrix2d S;
                                S(0, 0) = measCov(0, 0) + predCov(projector[0], projector[0]);
                                S(0, 1) = measCov(0, 1) + predCov(projector[0], projector[1]);
                                S(1, 0) = measCov(1, 0) + predCov(projector[1], projector[0]);
                                S(1, 1) = measCov(1, 1) + predCov(projector[1], projector[1]);
                                if (S.allFinite() && S(0, 0) > 0.0 && S(1, 1) > 0.0) {
                                    innovPull0 = innov0 / std::sqrt(S(0, 0));
                                    innovPull1 = innov1 / std::sqrt(S(1, 1));
                                    const Eigen::Vector2d r(innov0, innov1);
                                    const auto ldlt = S.ldlt();
                                    if (ldlt.info() == Eigen::Success) innovChi2 = r.dot(ldlt.solve(r));
                                }
                            }
                        } catch (...) {}
                    }
                }

""" + innov_anchor
cc = replace_once(cc, innov_anchor, innov_code, "predicted innovation")

cc = replace_once(
    cc,
    "                    ++m_nSensorMapExact;\n",
    "                    ++out.nMapExact;\n",
    "exact map counter",
)
cc = replace_once(
    cc,
    "                        ++m_nSensorMapFallback;\n",
    "                        ++out.nMapFallback;\n",
    "fallback map counter",
)
cc = replace_once(
    cc,
    "                        ++m_nSensorMapFailed;\n                    }\n                } else {\n                    ++m_nSensorMapFailed;\n                }",
    "                        ++out.nMapFailedPhysics;\n                    }\n                } else {\n                    ++out.nMapUnmappedNonPhysics;\n                }",
    "failed map counters",
)

push_anchor = """                out.state_track_index.push_back(actsTrackIndex);
                out.state_index.push_back(stateIndex++);
"""
push_new = """                out.state_track_index.push_back(actsTrackIndex);
                out.state_parent_seed_index.push_back(parentSeedIndex);
                out.state_parent_track_index.push_back(parentTrackIndex);
                out.state_parent_track_collectionID.push_back(parentTrackCollectionID);
                out.state_parent_identity_valid.push_back(parentIdentityValid);
                out.state_estimate_kind.push_back(estimateKind);
                out.state_index.push_back(stateIndex++);
"""
cc = replace_once(cc, push_anchor, push_new, "state identity push")
cc = replace_once(
    cc,
    "                out.state_resid_loc1.push_back(resid1);\n                out.x_on_plane.push_back(globalX);",
    "                out.state_resid_loc1.push_back(resid1);\n"
    "                out.state_meas_dim.push_back(measDim);\n"
    "                out.state_proj_index0.push_back(proj0); out.state_proj_index1.push_back(proj1);\n"
    "                out.state_meas_cov00.push_back(measCov00); out.state_meas_cov01.push_back(measCov01);\n"
    "                out.state_meas_cov11.push_back(measCov11);\n"
    "                out.state_pred_loc0.push_back(predLoc0); out.state_pred_loc1.push_back(predLoc1);\n"
    "                out.state_pred_theta.push_back(predTheta); out.state_pred_phi.push_back(predPhi);\n"
    "                out.state_pred_qOverP.push_back(predQOverP); out.state_pred_time.push_back(predTime);\n"
    "                for (std::size_t k = 0; k < predCovUpper.size(); ++k) out.state_pred_cov_upper[k].push_back(predCovUpper[k]);\n"
    "                out.state_innov0.push_back(innov0); out.state_innov1.push_back(innov1);\n"
    "                out.state_innov_pull0.push_back(innovPull0); out.state_innov_pull1.push_back(innovPull1);\n"
    "                out.state_innov_chi2.push_back(innovChi2);\n"
    "                out.x_on_plane.push_back(globalX);",
    "innovation push",
)
cc = replace_once(
    cc,
    "                out.state_pdg.push_back(\n                    (actsTrackIndex >= 0 && static_cast<std::size_t>(actsTrackIndex) < out.pdg.size())\n                        ? out.pdg[static_cast<std::size_t>(actsTrackIndex)]\n                        : 0);",
    "                int statePdg = 0;\n"
    "                if (parentIdentityValid) {\n"
    "                    if (const auto found = pdgByTrackObject.find({parentTrackCollectionID, parentTrackIndex});\n"
    "                        found != pdgByTrackObject.end()) statePdg = found->second;\n"
    "                }\n"
    "                out.state_pdg.push_back(statePdg);",
    "state pdg stable",
)

cc = replace_once(
    cc,
    "    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks, ckfActsTrackStates);\n\n    m_tree->Fill();",
    "    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks, ckfActsTrackStates);\n"
    "    m_selPrimaryHasTruthMatchedTrack = m_ckf.truthMatched.index >= 0 ? 1 : 0;\n"
    "    m_nSensorMapExact = m_ts.nMapExact + m_ckf.nMapExact;\n"
    "    m_nSensorMapFallback = m_ts.nMapFallback + m_ckf.nMapFallback;\n"
    "    m_nSensorMapFailed = m_ts.nMapFailedPhysics + m_ckf.nMapFailedPhysics;\n\n"
    "    m_tree->Fill();",
    "final stage/mapping totals",
)

# Clear/bind new fields.
cc = replace_once(
    cc,
    "    assocMcIndex = -1;\n    pdg = 0;",
    "    assocMcIndex = -1;\n    assocMcCollectionID = 0;\n    objectIndex = -1; objectCollectionID = 0; seedIndex = -1; identityValid = 0;\n    pdg = 0;",
    "BestSel reset",
)
cc = replace_once(
    cc,
    "    index.clear(); charge.clear(); type.clear(); pdg.clear();\n",
    "    index.clear(); charge.clear(); type.clear(); pdg.clear();\n"
    "    object_index.clear(); object_collectionID.clear(); seed_index.clear(); identity_valid.clear();\n"
    "    for (auto& v : cov_upper) v.clear();\n",
    "track clear identity",
)
cc = replace_once(
    cc,
    "    state_track_index.clear(); state_index.clear(); state_acts_index.clear();\n",
    "    state_track_index.clear(); state_index.clear(); state_acts_index.clear();\n"
    "    state_parent_seed_index.clear(); state_parent_track_index.clear();\n"
    "    state_parent_track_collectionID.clear(); state_parent_identity_valid.clear();\n"
    "    state_estimate_kind.clear();\n",
    "state clear identity",
)
cc = replace_once(
    cc,
    "    state_resid_loc0.clear(); state_resid_loc1.clear();\n",
    "    state_resid_loc0.clear(); state_resid_loc1.clear();\n"
    "    state_meas_dim.clear(); state_proj_index0.clear(); state_proj_index1.clear();\n"
    "    state_meas_cov00.clear(); state_meas_cov01.clear(); state_meas_cov11.clear();\n"
    "    state_pred_loc0.clear(); state_pred_loc1.clear(); state_pred_theta.clear(); state_pred_phi.clear();\n"
    "    state_pred_qOverP.clear(); state_pred_time.clear();\n"
    "    for (auto& v : state_pred_cov_upper) v.clear();\n"
    "    state_innov0.clear(); state_innov1.clear(); state_innov_pull0.clear(); state_innov_pull1.clear();\n"
    "    state_innov_chi2.clear();\n",
    "innovation clear",
)
cc = replace_once(
    cc,
    "    hasTrack = 0;\n}",
    "    nMapExact = 0; nMapFallback = 0; nMapFailedPhysics = 0; nMapUnmappedNonPhysics = 0;\n"
    "    hasTrack = 0;\n}",
    "map clear",
)

cc = replace_once(
    cc,
    "    br(prefix + \"assoc_mcIndex\", &b.assocMcIndex);\n    br(prefix + \"assoc_weight\", &b.assocWeight);",
    "    br(prefix + \"assoc_mcIndex\", &b.assocMcIndex);\n"
    "    br(prefix + \"assoc_mcCollectionID\", &b.assocMcCollectionID);\n"
    "    br(prefix + \"object_index\", &b.objectIndex);\n"
    "    br(prefix + \"object_collectionID\", &b.objectCollectionID);\n"
    "    br(prefix + \"seed_index\", &b.seedIndex);\n"
    "    br(prefix + \"identity_valid\", &b.identityValid);\n"
    "    br(prefix + \"assoc_weight\", &b.assocWeight);",
    "best bind identity",
)
cc = replace_once(
    cc,
    "    br(trkPrefix + \"index\", &c.index);\n    br(trkPrefix + \"type\", &c.type);",
    "    br(trkPrefix + \"index\", &c.index);\n"
    "    br(trkPrefix + \"object_index\", &c.object_index);\n"
    "    br(trkPrefix + \"object_collectionID\", &c.object_collectionID);\n"
    "    br(trkPrefix + \"seed_index\", &c.seed_index);\n"
    "    br(trkPrefix + \"identity_valid\", &c.identity_valid);\n"
    "    {\n"
    "        static const std::array<std::pair<int,int>,21> ij{{\n"
    "            {0,0},{0,1},{0,2},{0,3},{0,4},{0,5},{1,1},{1,2},{1,3},{1,4},{1,5},\n"
    "            {2,2},{2,3},{2,4},{2,5},{3,3},{3,4},{3,5},{4,4},{4,5},{5,5}}};\n"
    "        for (std::size_t k = 0; k < ij.size(); ++k)\n"
    "            br(trkPrefix + \"cov_\" + std::to_string(ij[k].first) + std::to_string(ij[k].second), &c.cov_upper[k]);\n"
    "    }\n"
    "    br(trkPrefix + \"type\", &c.type);",
    "track bind identity/cov",
)
cc = replace_once(
    cc,
    "    br(trkPrefix + \"state_track_index\", &c.state_track_index);\n    br(trkPrefix + \"state_index\", &c.state_index);",
    "    br(trkPrefix + \"state_track_index\", &c.state_track_index);\n"
    "    br(trkPrefix + \"state_parent_seed_index\", &c.state_parent_seed_index);\n"
    "    br(trkPrefix + \"state_parent_track_index\", &c.state_parent_track_index);\n"
    "    br(trkPrefix + \"state_parent_track_collectionID\", &c.state_parent_track_collectionID);\n"
    "    br(trkPrefix + \"state_parent_identity_valid\", &c.state_parent_identity_valid);\n"
    "    br(trkPrefix + \"state_estimate_kind\", &c.state_estimate_kind);\n"
    "    br(trkPrefix + \"state_index\", &c.state_index);",
    "state bind identity",
)
cc = replace_once(
    cc,
    "    br(trkPrefix + \"state_resid_loc1\", &c.state_resid_loc1);\n    br(trkPrefix + \"x_on_plane\", &c.x_on_plane);",
    "    br(trkPrefix + \"state_resid_loc1\", &c.state_resid_loc1);\n"
    "    br(trkPrefix + \"state_meas_dim\", &c.state_meas_dim);\n"
    "    br(trkPrefix + \"state_proj_index0\", &c.state_proj_index0); br(trkPrefix + \"state_proj_index1\", &c.state_proj_index1);\n"
    "    br(trkPrefix + \"state_meas_cov00\", &c.state_meas_cov00); br(trkPrefix + \"state_meas_cov01\", &c.state_meas_cov01);\n"
    "    br(trkPrefix + \"state_meas_cov11\", &c.state_meas_cov11);\n"
    "    br(trkPrefix + \"state_pred_loc0\", &c.state_pred_loc0); br(trkPrefix + \"state_pred_loc1\", &c.state_pred_loc1);\n"
    "    br(trkPrefix + \"state_pred_theta\", &c.state_pred_theta); br(trkPrefix + \"state_pred_phi\", &c.state_pred_phi);\n"
    "    br(trkPrefix + \"state_pred_qOverP\", &c.state_pred_qOverP); br(trkPrefix + \"state_pred_time\", &c.state_pred_time);\n"
    "    {\n"
    "        static const std::array<std::pair<int,int>,21> ij{{\n"
    "            {0,0},{0,1},{0,2},{0,3},{0,4},{0,5},{1,1},{1,2},{1,3},{1,4},{1,5},\n"
    "            {2,2},{2,3},{2,4},{2,5},{3,3},{3,4},{3,5},{4,4},{4,5},{5,5}}};\n"
    "        for (std::size_t k = 0; k < ij.size(); ++k)\n"
    "            br(trkPrefix + \"state_pred_cov_\" + std::to_string(ij[k].first) + std::to_string(ij[k].second), &c.state_pred_cov_upper[k]);\n"
    "    }\n"
    "    br(trkPrefix + \"state_innov0\", &c.state_innov0); br(trkPrefix + \"state_innov1\", &c.state_innov1);\n"
    "    br(trkPrefix + \"state_innov_pull0\", &c.state_innov_pull0); br(trkPrefix + \"state_innov_pull1\", &c.state_innov_pull1);\n"
    "    br(trkPrefix + \"state_innov_chi2\", &c.state_innov_chi2);\n"
    "    br(trkPrefix + \"x_on_plane\", &c.x_on_plane);",
    "innovation bind",
)
cc = replace_once(
    cc,
    "    br(trkPrefix + \"has_track\", &c.hasTrack);\n",
    "    br(trkPrefix + \"n_sensor_map_exact\", &c.nMapExact);\n"
    "    br(trkPrefix + \"n_sensor_map_fallback\", &c.nMapFallback);\n"
    "    br(trkPrefix + \"n_sensor_map_failed_physics\", &c.nMapFailedPhysics);\n"
    "    br(trkPrefix + \"n_sensor_map_unmapped_nonphysics\", &c.nMapUnmappedNonPhysics);\n"
    "    br(trkPrefix + \"has_track\", &c.hasTrack);\n",
    "map bind",
)
ccpath.write_text(cc)

# Add pure helper tests for the schema-changing semantics.
tpath = Path("tests/test_helpers.cc")
t = tpath.read_text()
t = replace_once(
    t,
    "    CHECK(closeTo(wrapPi(0.1), 0.1));\n",
    "    CHECK(closeTo(wrapPi(0.1), 0.1));\n\n"
    "    CHECK(kSchemaVersion == 3);\n"
    "    CHECK(preferredEstimateKind(true, false, false) == kEstimatePredicted);\n"
    "    CHECK(preferredEstimateKind(true, true, false) == kEstimateFiltered);\n"
    "    CHECK(preferredEstimateKind(true, true, true) == kEstimateSmoothed);\n"
    "    CHECK(preferredEstimateKind(false, false, false) == kEstimateNone);\n"
    "    CHECK(closeTo(nativeTimeToNs(299.792458, 299.792458), 1.0));\n"
    "    CHECK(std::isnan(nativeTimeToNs(1.0, 0.0)));\n",
    "schema3 helper tests",
)
tpath.write_text(t)
