# Inverse optical flow
Python library for computing inverse optical flow.

Sánchez, J., Salgado, A., Monzón, N. (2013). An Efficient Algorithm for Estimating the Inverse Optical Flow. In: Sanches, J.M., Micó, L., Cardoso, J.S. (eds) Pattern Recognition and Image Analysis. IbPRIA 2013. Lecture Notes in Computer Science, vol 7887. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-38628-2_46

[./inverse_flow](`./inverse_flow`) folder contains original source code from article attachments.


## Installation
```shell
pip install git+https://github.com/ivanstepanovftw/inverse_optical_flow.git@main
```

## Usage
```python
import numpy as np
from inverse_optical_flow import inverse_flow_max, inverse_flow_avg

# Shape of flow is (height, width, 2)
forward_flow = np.array([
    [[0, 0], [0, 2], [0, 0]],
    [[0, 0], [1, 1], [0, 0]],
    [[0, 0], [0, 0], [0, 0]],
], dtype=np.float32)

# Compute backward flow using max method
backward_flow, no_disocclusion_mask = inverse_flow_max(forward_flow)

# See results
assert np.allclose(backward_flow, np.array([
    [[0, 0], [0,  0], [ 0,  0]],
    [[0, 0], [0,  0], [ 0,  0]],
    [[0, 0], [0, -2], [-1, -1]],
]))
assert np.allclose(no_disocclusion_mask, np.array([
    [255,   0, 255],
    [255,   0, 255],
    [255, 255, 255]
]))
```

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without any
additional terms or conditions.
