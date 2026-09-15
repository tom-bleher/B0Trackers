from pathlib import Path

path = Path("B0Trackers.cc")
text = path.read_text()

old_tp = '''            auto tp = (trajectory->trackParameters_size() > 0)
                ? trajectory->getTrackParameters(0)
                : edm4eic::TrackParameters::makeEmpty();
            if (!tp.isAvailable() && trajIndex < tracks.size() && tracks[trajIndex] != nullptr) {
                tp = *tracks[trajIndex];
            }
            if (!tp.isAvailable()) continue;
'''
new_tp = '''            auto tp = (trajectory->trackParameters_size() > 0)
                ? trajectory->getTrackParameters(0)
                : edm4eic::TrackParameters::makeEmpty();
            // Schema 3 never guesses trajectory identity from parallel collection
            // positions. If the trajectory does not carry fitted parameters, leave
            // this object unresolved rather than borrowing tracks[trajIndex].
            if (!tp.isAvailable()) continue;
'''
if old_tp not in text:
    raise SystemExit("TrackParameters positional fallback block not found")
text = text.replace(old_tp, new_tp, 1)

old_track = '''            if (stableEdmTrack != nullptr) {
                chi2 = stableEdmTrack->getChi2();
                ndf = static_cast<int>(stableEdmTrack->getNdf());
                pdg = stableEdmTrack->getPdg();
            } else if (trajIndex < edmTracks.size() && edmTracks[trajIndex] != nullptr) {
                chi2 = edmTracks[trajIndex]->getChi2();
                ndf = static_cast<int>(edmTracks[trajIndex]->getNdf());
                pdg = edmTracks[trajIndex]->getPdg();
                const auto objectId = edmTracks[trajIndex]->id();
                objectIndex = objectId.index;
                objectCollectionID = objectId.collectionID;
            }
'''
new_track = '''            if (stableEdmTrack != nullptr) {
                chi2 = stableEdmTrack->getChi2();
                ndf = static_cast<int>(stableEdmTrack->getNdf());
                pdg = stableEdmTrack->getPdg();
            }
'''
if old_track not in text:
    raise SystemExit("EDM track positional fallback block not found")
text = text.replace(old_track, new_track, 1)

path.write_text(text)
