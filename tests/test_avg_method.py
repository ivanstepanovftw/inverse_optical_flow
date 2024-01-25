from math import nan

import numpy as np
import inverse_optical_flow

# Shape of the flow is (2, height, width).
# The first channel is horizontal, x-axis flow. The second channel is vertical, y-axis flow.
forward_flow = np.array([
    [[0, 0, 0],
     [0, 1, 0],
     [0, 0, 0]],

    [[0, 2, 0],
     [0, 1, 0],
     [0, 0, 0]],
], dtype=np.float32)

# Compute backward flow using avg method
backward_flow, disocclusion_mask = inverse_optical_flow.avg_method(forward_flow)

# See results
assert np.allclose(backward_flow, np.array([
    [[0,  0,   0],
     [0,  0,   0],
     [0,  0, nan]],

    [[0,  0,   0],
     [0,  0,   0],
     [0, -2, nan]]
]), equal_nan=True), backward_flow

assert np.allclose(disocclusion_mask, np.array([
    [1, 1, 1],
    [1, 1, 1],
    [1, 0, 0]
])), disocclusion_mask
