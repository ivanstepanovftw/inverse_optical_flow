from glob import glob
# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.2"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

setup(
    name="inverse_optical_flow",
    version=__version__,
    description="Estimate inverse optical flow and disocclusion mask from optical flow.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Ivan Stepanov",
    author_email="ivanstepanovftw@gmail.com",
    url="https://github.com/ivanstepanovftw/inverse_optical_flow",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Libraries",
        "License :: OSI Approved :: MIT License",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="optical flow, inverse optical flow, computer vision, image processing, video processing",
    ext_modules=[
        Pybind11Extension(
            "inverse_optical_flow",
            sorted(glob("src/*.cpp")),
            # Example: passing in the version to the compiled code
            define_macros=[('VERSION_INFO', __version__)],
        ),
    ],
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    install_requires=[
        'numpy',
    ],
    python_requires='>=3.7, <4',
    package_data={
        'inverse_optical_flow': ['LICENSE-APACHE-2.0', 'LICENSE-MIT', 'README.md'],
    },
)
