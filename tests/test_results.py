import json

import numpy as np

from pstokes_fem.metrics import ERROR_KEYS
from pstokes_fem.results import load_result, save_result


def test_npz_roundtrip_without_pickle(tmp_path):
    series = {key: [1.0, 2.0] for key in ERROR_KEYS}
    metadata = {
        "level": 0,
        "time_step": 0.025,
        "parameters": {"power_law": 1.5},
    }
    path = save_result(tmp_path / "result.npz", series, metadata)
    loaded, loaded_metadata = load_result(path)
    assert loaded_metadata == metadata
    for key in ERROR_KEYS:
        np.testing.assert_allclose(loaded[key], series[key])

