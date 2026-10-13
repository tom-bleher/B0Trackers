# B0Trackers schema 7: truth-seeded unfiltered candidates

Schema 7 adds the truth-seeded chain to the `ckfdiag_unf_*` rows
(`B0TrackerCKFTruthSeededActsTrackStatesUnfiltered` /
`B0TrackerCKFTruthSeededActsTracksUnfiltered`, read only with
`write_track_states=1`).

## Branches

- `ckfdiag_unf_chain`: 0 = stub-seeded, 1 = truth-seeded.
- `ckfdiag_unf_mcindex` / `ckfdiag_unf_mccol`: per-row MC association. Stub
  rows carry the selected primary (-1 when none). Truth-seeded rows carry the
  particle from `B0TrackerCKFTruthSeededTrackUnfilteredAssociations` (-1 when
  the candidate has no unfiltered EDM track).
- `has_ts_unfiltered_assocs`, `has_ts_acts_states_unfiltered`,
  `has_ts_acts_tracks_unfiltered`: collection presence flags.

All other `ckfdiag_unf_*` branches are shared. For truth-seeded rows the
truth residuals anchor on the row's own particle (`tdist_prim` =
distance to own hits, `tdist_other` = nearest other particle).

## Joining truth-seeded rows to primaries

`ckfdiag_unf_seedpos` is the truth-seed position for chain 1. The
per-truth-seed `truth_seed_p/theta/phi/charge` branches (parallel to the
truth-seed input order) identify the particle kinematically when the
EDM association is absent: match charge exactly and theta/phi closely
against the selected primary (p within the truth-seed momentum smear).
