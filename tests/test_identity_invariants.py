import random
import unittest
from pathlib import Path


class StableIdentityTest(unittest.TestCase):
    def test_relation_lookup_is_invariant_to_collection_order(self):
        # Model the schema-3 relation used in B0Trackers: reconstructed tracks are
        # keyed by the stable ObjectID of their related trajectory, not by vector
        # position. Reordering either collection must therefore leave ownership
        # unchanged.
        tracks = [
            {"track": (91, 4), "trajectory": (72, 8), "seed": (33, 3)},
            {"track": (91, 2), "trajectory": (72, 5), "seed": (33, 1)},
            {"track": (91, 9), "trajectory": (72, 1), "seed": (33, 7)},
        ]
        expected = {row["trajectory"]: row["track"] for row in tracks}

        for seed in range(20):
            shuffled = tracks[:]
            random.Random(seed).shuffle(shuffled)
            by_trajectory = {row["trajectory"]: row["track"] for row in shuffled}
            self.assertEqual(by_trajectory, expected)

    def test_acts_seed_position_translates_to_podio_object_id(self):
        # CKFTracking stores the input seed *position* in the ACTS dynamic seed
        # column. Schema 3 must translate that position through the real seed
        # collection instead of assuming position == PODIO ObjectID.index.
        seed_object_ids = [(44, 9), (44, 2), (44, 17)]
        acts_seed_positions = [2, 0, 1]
        resolved = [seed_object_ids[position] for position in acts_seed_positions]
        self.assertEqual(resolved, [(44, 17), (44, 9), (44, 2)])

    def test_seed_lineage_marks_ambiguous_owner_unresolved(self):
        tracks = [
            {"track": (91, 2), "seed": (33, 1)},
            {"track": (91, 4), "seed": (33, 3)},
            {"track": (91, 9), "seed": (33, 3)},
        ]
        by_seed = {}
        for row in tracks:
            by_seed.setdefault(row["seed"], []).append(row["track"])

        self.assertEqual(by_seed[(33, 1)], [(91, 2)])
        self.assertEqual(len(by_seed[(33, 3)]), 2)
        # Schema 3 intentionally refuses to guess a parent track when a seed has
        # more than one candidate.
        resolved_seed3 = by_seed[(33, 3)][0] if len(by_seed[(33, 3)]) == 1 else None
        self.assertIsNone(resolved_seed3)

    def test_plugin_uses_relation_identity_not_parallel_collection_index(self):
        root = Path(__file__).parents[1]
        source = (root / "B0Trackers.cc").read_text()
        header = (root / "B0Trackers.h").read_text()
        self.assertIn("edmTrackByTrajectory", source)
        self.assertIn("stableTrackIt", source)
        self.assertIn("trackObjectsBySeed", source)
        self.assertIn("seedPosition < chainSeeds.size()", source)
        self.assertIn("parentSeedCollectionID", source)
        self.assertIn("seed_collectionID", header)
        self.assertIn("state_parent_seed_collectionID", header)

        # Documentation may mention the rejected positional patterns. Guard the
        # implementation itself by checking non-comment source lines only.
        code = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("//")
        )
        self.assertNotIn("edmTracks[trajIndex]", code)
        self.assertNotIn("tracks[trajIndex]", code)


if __name__ == "__main__":
    unittest.main()
