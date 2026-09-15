from pathlib import Path

h = Path("B0Trackers.h")
text = h.read_text()
text = text.replace(
    "    int m_primaryPdg = 2212;\n    int m_primaryStatus = 1;",
    "    int m_primaryPdg = 2212;\n    int m_primaryStatus = 1;\n    int m_minMeasurementStations = 3;",
    1,
)
text = text.replace(
    "    int m_nStationsPrimary = 0; // selected primary only",
    "    int m_nStationsPrimary = 0; // selected primary only\n    int m_nSelectedPrimaryMeasurements = -1;\n    int m_nMeasurementStationsSelectedPrimary = -1;\n    int m_selPrimaryMeasurementReconstructable = -1;",
    1,
)
h.write_text(text)

p = Path("B0Trackers.cc")
text = p.read_text()
text = text.replace(
    '''    app->SetDefaultParameter("B0Trackers:primary_status", m_primaryStatus,
                             "Generator status used to tag the selected primary particle");
''',
    '''    app->SetDefaultParameter("B0Trackers:primary_status", m_primaryStatus,
                             "Generator status used to tag the selected primary particle");
    app->SetDefaultParameter("B0Trackers:min_measurement_stations", m_minMeasurementStations,
                             "Minimum distinct selected-primary measurement stations for reconstructability");
''',
    1,
)
text = text.replace(
    '''    m_tree->Branch("n_stations_primary", &m_nStationsPrimary);
''',
    '''    m_tree->Branch("n_stations_primary", &m_nStationsPrimary);
    m_tree->Branch("min_measurement_stations_required", &m_minMeasurementStations);
    m_tree->Branch("n_measurements_selected_primary", &m_nSelectedPrimaryMeasurements);
    m_tree->Branch("n_measurement_stations_selected_primary",
                   &m_nMeasurementStationsSelectedPrimary);
    m_tree->Branch("sel_primary_measurement_reconstructable",
                   &m_selPrimaryMeasurementReconstructable);
''',
    1,
)
text = text.replace(
    '''    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;
    {
        // (cellID, MC ObjectID) -> summed deposit.
        std::map<std::pair<std::uint64_t, std::pair<std::uint32_t, int>>, double> cellParticleEDep;
''',
    '''    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;
    // Full per-cell truth composition is retained for measurement-level
    // reconstructability; SimLink below is only the dominant-cell summary.
    std::map<std::uint64_t, std::map<std::pair<std::uint32_t, int>, double>> cellTruthEDep;
    {
        // (cellID, MC ObjectID) -> summed deposit.
        std::map<std::pair<std::uint64_t, std::pair<std::uint32_t, int>>, double> cellParticleEDep;
''',
    1,
)
text = text.replace(
    '''            cellParticleEDep[{cid, {id.collectionID, id.index}}] += sim.getEDep();
            auto& link = simLinkByCell[cid];
''',
    '''            cellParticleEDep[{cid, {id.collectionID, id.index}}] += sim.getEDep();
            cellTruthEDep[cid][{id.collectionID, id.index}] += sim.getEDep();
            auto& link = simLinkByCell[cid];
''',
    1,
)
anchor = '''    const auto isSelPrimary = [this](std::uint32_t collectionID, int index) -> int {
        return (m_selPrimaryMcIndex >= 0 &&
                collectionID == m_selPrimaryMcCollectionID &&
                index == m_selPrimaryMcIndex) ? 1 : 0;
    };

'''
insert = anchor + '''    // A truth crossing is only geometrical acceptance. Reconstructability also
    // requires digitized/reconstructed measurements attributable to the selected
    // primary in enough distinct physical B0 stations. Build the measurement
    // truth label from the summed truth composition of its constituent raw cells.
    m_nSelectedPrimaryMeasurements = hasRawAssocs ? 0 : -1;
    m_nMeasurementStationsSelectedPrimary = hasRawAssocs ? 0 : -1;
    m_selPrimaryMeasurementReconstructable = hasRawAssocs ? 0 : -1;
    if (hasRawAssocs && m_selPrimaryMcIndex >= 0) {
        std::set<int> selectedMeasurementStations;
        for (const auto* measurement : measurements) {
            if (measurement == nullptr) continue;
            std::map<std::pair<std::uint32_t, int>, double> measurementTruth;
            std::set<std::uint64_t> seenCells;
            std::set<int> measurementStations;
            for (const auto& hit : measurement->getHits()) {
                const auto raw = hit.getRawHit();
                if (!raw.isAvailable()) continue;
                const auto cid = static_cast<std::uint64_t>(raw.getCellID());
                if (!seenCells.insert(cid).second) continue;
                int plane = -1, module = -1, sensor = -1, side = -1;
                decodeIds(cid, plane, module, sensor, side);
                const int station = stationOf(plane);
                if (station > 0) measurementStations.insert(station);
                const auto truthIt = cellTruthEDep.find(cid);
                if (truthIt == cellTruthEDep.end()) continue;
                for (const auto& [mcId, edep] : truthIt->second) {
                    measurementTruth[mcId] += edep;
                }
            }
            if (measurementTruth.empty()) continue;
            const auto dominant = std::max_element(
                measurementTruth.begin(), measurementTruth.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            if (dominant == measurementTruth.end()) continue;
            if (dominant->first.first != m_selPrimaryMcCollectionID ||
                dominant->first.second != m_selPrimaryMcIndex) {
                continue;
            }
            ++m_nSelectedPrimaryMeasurements;
            selectedMeasurementStations.insert(measurementStations.begin(), measurementStations.end());
        }
        m_nMeasurementStationsSelectedPrimary =
            static_cast<int>(selectedMeasurementStations.size());
        m_selPrimaryMeasurementReconstructable =
            m_nMeasurementStationsSelectedPrimary >= m_minMeasurementStations ? 1 : 0;
    }

'''
if anchor not in text:
    raise SystemExit("selected-primary anchor not found")
text = text.replace(anchor, insert, 1)
text = text.replace(
    '''    m_nStationsPrimary = 0;
''',
    '''    m_nStationsPrimary = 0;
    m_nSelectedPrimaryMeasurements = hasRawAssocs ? 0 : -1;
    m_nMeasurementStationsSelectedPrimary = hasRawAssocs ? 0 : -1;
    m_selPrimaryMeasurementReconstructable = hasRawAssocs ? 0 : -1;
''',
    1,
)
p.write_text(text)
