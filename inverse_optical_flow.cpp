#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef WEIGHT_TH
#define WEIGHT_TH 0.25
#endif
#ifndef MOTION_TH
#define MOTION_TH 0.25
#endif

namespace py = pybind11;


auto inverse_flow_max(const py::array_t<float> & flow_array) -> std::pair<py::array_t<float>, py::array_t<uint8_t>> {
    auto flow = flow_array.unchecked<3>();
    const auto ny = flow.shape(0);
    const auto nx = flow.shape(1);
    if (flow.shape(2) != 2)
        throw std::runtime_error("Input flow must have shape (ny, nx, 2)");

    auto inverse_flow_array = py::array_t<float>({ny, nx, ssize_t(2)});
    auto no_disocclusion_mask_array = py::array_t<uint8_t>({ny, nx});
    inverse_flow_array[py::make_tuple(py::ellipsis())] = 0.f;
    no_disocclusion_mask_array[py::make_tuple(py::ellipsis())] = 0;
    auto flow_i = inverse_flow_array.mutable_unchecked<3>();
    auto no_disocclusion_mask = no_disocclusion_mask_array.mutable_unchecked<2>();

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            // warping the flow
            const auto yw = float(y) + flow(y, x, 1);
            const auto xw = float(x) + flow(y, x, 0);
            // integer part of the warped position
            auto yi = ssize_t(yw);
            auto xi = ssize_t(xw);
            // sign of the warped position
            const int sy = (yw < 0) ? -1 : 1;
            const int sx = (xw < 0) ? -1 : 1;
            // warped position
            ssize_t dy = yi + sy;
            ssize_t dx = xi + sx;
            // check that the warped position is inside the image
            yi = std::max(ssize_t(0), std::min(ny - 1, yi));
            xi = std::max(ssize_t(0), std::min(nx - 1, xi));
            dy = std::max(ssize_t(0), std::min(ny - 1, dy));
            dx = std::max(ssize_t(0), std::min(nx - 1, dx));
            // compute the four proportions
            const auto e1 = float(sx) * (xw - float(xi));
            const auto E1 = 1.0f - e1;
            const auto e2 = float(sy) * (yw - float(yi));
            const auto E2 = 1.0f - e2;
            // put in the four points the corresponding proportion
            const auto w1 = E1 * E2;
            const auto w2 = e1 * E2;
            const auto w3 = E1 * e2;
            const auto w4 = e1 * e2;
            // compute the four distances
            const float d  = std::pow(flow(y, x, 1), 2) + std::pow(flow(y, x, 0), 2);
            const float d1 = std::pow(flow_i(yi, xi, 1), 2) + std::pow(flow_i(yi, xi, 0), 2);
            const float d2 = std::pow(flow_i(yi, dx, 1), 2) + std::pow(flow_i(yi, dx, 0), 2);
            const float d3 = std::pow(flow_i(dy, xi, 1), 2) + std::pow(flow_i(dy, xi, 0), 2);
            const float d4 = std::pow(flow_i(dy, dx, 1), 2) + std::pow(flow_i(dy, dx, 0), 2);

            // check if the warped position is occluded
            if (w1 >= WEIGHT_TH && d >= d1) {
                flow_i(yi, xi, 0) = -flow(y, x, 0);
                flow_i(yi, xi, 1) = -flow(y, x, 1);
                no_disocclusion_mask(yi, xi) = 255;
            }

            if (w2 >= WEIGHT_TH && d >= d2) {
                flow_i(yi, dx, 0) = -flow(y, x, 0);
                flow_i(yi, dx, 1) = -flow(y, x, 1);
                no_disocclusion_mask(yi, dx) = 255;
            }

            if (w3 >= WEIGHT_TH && d >= d3) {
                flow_i(dy, xi, 0) = -flow(y, x, 0);
                flow_i(dy, xi, 1) = -flow(y, x, 1);
                no_disocclusion_mask(dy, xi) = 255;
            }

            if (w4 >= WEIGHT_TH && d >= d4) {
                flow_i(dy, dx, 0) = -flow(y, x, 0);
                flow_i(dy, dx, 1) = -flow(y, x, 1);
                no_disocclusion_mask(dy, dx) = 255;
            }
        }
    }

    return std::make_pair(inverse_flow_array, no_disocclusion_mask_array);
}


inline void select_motion(
    const float d,
    const float u,
    const float v,
    const float wght,
    float & d_,
    float & u_,
    float & v_,
    float & wght_,
    float & mask
) {
    if (wght >= WEIGHT_TH) {
        if (fabs(d - d_) <= MOTION_TH) {
            u_ += u * wght;
            v_ += v * wght;
            wght_ += wght;
            mask = 255;
        } else if (d >= d_) {
            //if it is an oclussion we retain the highest value
            d_ = d;
            u_ = u * wght;
            v_ = v * wght;
            wght_ = wght;
            mask = 255;
        }
    }
}

auto inverse_flow_avg(const py::array_t<float> & flow_array) -> std::pair<py::array_t<float>, py::array_t<uint8_t>> {
    auto flow = flow_array.unchecked<3>();
    const auto ny = flow.shape(0);
    const auto nx = flow.shape(1);
    if (flow.shape(2) != 2)
        throw std::runtime_error("Input flow must have shape (ny, nx, 2)");

    // Define the output arrays
    auto inverse_flow_array = py::array_t<float>({ny, nx, ssize_t(2)});
    auto no_disocclusion_mask_array = py::array_t<uint8_t>({ny, nx});
    inverse_flow_array[py::make_tuple(py::ellipsis())] = 0.f;
    no_disocclusion_mask_array[py::make_tuple(py::ellipsis())] = 0;

    // Temporary arrays
    auto avg_uv_array = py::array_t<float>({ny, nx, ssize_t(2)});
    auto wgt_array = py::array_t<float>({ny, nx});
    auto d_array = py::array_t<float>({ny, nx});
    avg_uv_array[py::make_tuple(py::ellipsis())] = 0.f;
    wgt_array[py::make_tuple(py::ellipsis())] = 0.f;
    d_array[py::make_tuple(py::ellipsis())] = 0.f;

    auto flow_i = inverse_flow_array.mutable_unchecked<3>();
    auto no_disocclusion_mask = no_disocclusion_mask_array.mutable_unchecked<2>();
    auto avg_uv_ = avg_uv_array.mutable_unchecked<3>();
    auto wgt_ = wgt_array.mutable_unchecked<2>();
    auto d_ = d_array.mutable_unchecked<2>();

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            // warping the flow
            const auto xw = float(x) + flow(y, x, 0);
            const auto yw = float(y) + flow(y, x, 1);
            // integer part of the warped position
            auto xi = ssize_t(xw);
            auto yi = ssize_t(yw);
            // sign of the warped position
            const int sx = (xw < 0) ? -1 : 1;
            const int sy = (yw < 0) ? -1 : 1;
            // warped position
            ssize_t dx = xi + sx;
            ssize_t dy = yi + sy;
            // check that the warped position is inside the image
            xi = std::max(ssize_t(0), std::min(nx - 1, xi));
            yi = std::max(ssize_t(0), std::min(ny - 1, yi));
            dx = std::max(ssize_t(0), std::min(nx - 1, dx));
            dy = std::max(ssize_t(0), std::min(ny - 1, dy));
            // compute the four proportions
            const auto e1 = float(sx) * (xw - float(xi));
            const auto E1 = 1.0f - e1;
            const auto e2 = float(sy) * (yw - float(yi));
            const auto E2 = 1.0f - e2;
            // put in the four points the corresponding proportion
            const auto w1 = E1 * E2;
            const auto w2 = e1 * E2;
            const auto w3 = E1 * e2;
            const auto w4 = e1 * e2;
            // compute the four distances
            const auto d  = float(std::pow(flow(y, x, 0), 2) + std::pow(flow(y, x, 1), 2));
            // select motion
            if (d >= WEIGHT_TH) {
                if (fabs(d - d_(yi, xi)) <= MOTION_TH) {
                    avg_uv_(yi, xi, 0) += flow(y, x, 0) * w1;
                    avg_uv_(yi, xi, 1) += flow(y, x, 1) * w1;
                    wgt_(yi, xi) += w1;
                    no_disocclusion_mask(yi, xi) = 255;
                } else if (d >= d_(yi, xi)) {
                    //if it is an oclussion we retain the highest value
                    d_(yi, xi) = d;
                    avg_uv_(yi, xi, 0) = flow(y, x, 0) * w1;
                    avg_uv_(yi, xi, 1) = flow(y, x, 1) * w1;
                    wgt_(yi, xi) = w1;
                    no_disocclusion_mask(yi, xi) = 255;
                }
            }
            if (d >= WEIGHT_TH) {
                if (fabs(d - d_(yi, dx)) <= MOTION_TH) {
                    avg_uv_(yi, dx, 0) += flow(y, x, 0) * w2;
                    avg_uv_(yi, dx, 1) += flow(y, x, 1) * w2;
                    wgt_(yi, dx) += w2;
                    no_disocclusion_mask(yi, dx) = 255;
                } else if (d >= d_(yi, dx)) {
                    //if it is an oclussion we retain the highest value
                    d_(yi, dx) = d;
                    avg_uv_(yi, dx, 0) = flow(y, x, 0) * w2;
                    avg_uv_(yi, dx, 1) = flow(y, x, 1) * w2;
                    wgt_(yi, dx) = w2;
                    no_disocclusion_mask(yi, dx) = 255;
                }
            }
            if (d >= WEIGHT_TH) {
                if (fabs(d - d_(dy, xi)) <= MOTION_TH) {
                    avg_uv_(dy, xi, 0) += flow(y, x, 0) * w3;
                    avg_uv_(dy, xi, 1) += flow(y, x, 1) * w3;
                    wgt_(dy, xi) += w3;
                    no_disocclusion_mask(dy, xi) = 255;
                } else if (d >= d_(dy, xi)) {
                    //if it is an oclussion we retain the highest value
                    d_(dy, xi) = d;
                    avg_uv_(dy, xi, 0) = flow(y, x, 0) * w3;
                    avg_uv_(dy, xi, 1) = flow(y, x, 1) * w3;
                    wgt_(dy, xi) = w3;
                    no_disocclusion_mask(dy, xi) = 255;
                }
            }
            if (d >= WEIGHT_TH) {
                if (fabs(d - d_(dy, dx)) <= MOTION_TH) {
                    avg_uv_(dy, dx, 0) += flow(y, x, 0) * w4;
                    avg_uv_(dy, dx, 1) += flow(y, x, 1) * w4;
                    wgt_(dy, dx) += w4;
                    no_disocclusion_mask(dy, dx) = 255;
                } else if (d >= d_(dy, dx)) {
                    //if it is an oclussion we retain the highest value
                    d_(dy, dx) = d;
                    avg_uv_(dy, dx, 0) = flow(y, x, 0) * w4;
                    avg_uv_(dy, dx, 1) = flow(y, x, 1) * w4;
                    wgt_(dy, dx) = w4;
                    no_disocclusion_mask(dy, dx) = 255;
                }
            }
        }
    }

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            if (no_disocclusion_mask(y, x) == 255) {
                flow_i(y, x, 0) = -avg_uv_(y, x, 0) / wgt_(y, x);
                flow_i(y, x, 1) = -avg_uv_(y, x, 1) / wgt_(y, x);
            }
        }
    }

    return std::make_pair(inverse_flow_array, no_disocclusion_mask_array);
}

PYBIND11_MODULE(inverse_optical_flow, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: inverse_optical_flow

        .. autosummary::
           :toctree: _generate

           inverse_flow_max
           inverse_flow_avg
    )pbdoc";
    m.def("inverse_flow_max", &inverse_flow_max, py::arg().noconvert(), "Inverse flow using max distance");
    m.def("inverse_flow_avg", &inverse_flow_avg, py::arg().noconvert(), "Inverse flow averaging closest points");
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
