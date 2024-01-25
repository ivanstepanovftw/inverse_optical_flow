# Inverse Optical Flow

[![PyPI version](https://badge.fury.io/py/inverse-optical-flow.svg)](https://badge.fury.io/py/inverse-optical-flow)

A Python library for estimating inverse optical flow and the disocclusion mask from optical flow.

Unofficial implementation of "An Efficient Algorithm for Estimating the Inverse Optical Flow" ([public full-text](https://www.researchgate.net/publication/258547558_An_Efficient_Algorithm_for_Estimating_the_Inverse_Optical_Flow)).

[`/inverse_flow`](/inverse_flow) folder contains original source code from the paper. This library does not have implemented filling disocclusion. For this purpose, you can use original implementation from the paper at [`/inverse_flow/fill_disocclusions.h`](inverse_flow/fill_disocclusions.h).

## Glossary

- **Optical flow** - the pattern of apparent motion of objects, surfaces, or edges in a visual scene, caused by the relative motion between an observer and the scene.

    Optical flow is used in computer vision to detect and describe the motion of objects and surfaces in sequential images. It represents the displacement of points between two consecutive frames, caused by the movement of the object or the camera.

- **Forward flow** - optical flow computed from the $I_1$ to $I_2$ image.

    Forward flow represents the movement of pixels from the $I_1$ to the $I_2$ image.

- **Backward flow** - optical flow computed from the $I_2$ to the $I_1$ image.

    Backward flow traces the origins of pixels from image $I_2$ back to their positions in image $I_1$.

- **Inverse optical flow** - a technique or process that reverses the direction of optical flow.

    While traditional optical flow estimates the movement of pixels from a first image to a second image (forward flow), inverse optical flow inverts this process. It aims to estimate the movement from the second image back to the first (backward flow), essentially reversing the flow direction.

- **Disocclusion mask** - a mask identifying regions of an image that are visible in one frame but not in the next due to occlusion.

    Disocclusion masks are crucial in optical flow for identifying areas where objects have moved in or out of the frame, thereby revealing or hiding parts of the scene.

![https://imgur.com/a/DdKKDVr](https://i.imgur.com/4zX85U3.png)
Middlebury test sequences: The first column shows the source image; the second, the ground truth; the third, the backward flow (disocclusions in white color); and the fourth, the inverse flow.

## Installation

### Prerequisites

- C++ compiler
- Python 3.7+
- pip

### Install

```shell
pip install -U inverse_optical_flow
```

## Usage

```python
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

# Compute backward flow using max method
backward_flow, disocclusion_mask = inverse_optical_flow.max_method(forward_flow)

# See results
assert np.allclose(backward_flow, np.array([
    [[0,  0,  0],
     [0,  0,  0],
     [0,  0, -1]],

    [[0,  0,  0],
     [0,  0,  0],
     [0, -2, -1]]
])), backward_flow

assert np.allclose(disocclusion_mask, np.array([
    [0, 1, 0],
    [0, 1, 0],
    [0, 0, 0]
])), disocclusion_mask
```

## Alternatives

https://github.com/sniklaus/softmax-splatting

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any
additional terms or conditions.

## Acknowledgements

Sánchez, J., Salgado, A., Monzón, N. (2013). An Efficient Algorithm for Estimating the Inverse Optical Flow. In: Sanches, J.M., Micó, L., Cardoso, J.S. (eds) Pattern Recognition and Image Analysis. IbPRIA 2013. Lecture Notes in Computer Science, vol 7887. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-38628-2_46
