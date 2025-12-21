# Latitude/Longitude to UTM/UPS

This repository provides a set of C++ routines to convert coordinates between latitude/longitude and the UTM/UPS grids.

__Prior to 2025-12-21__ these were based on the series expansion provided by the US Defense Mapping Agency Technical Manual 8358.2.

__After 2025-12-21__ the UTM routines use an implementation based on the expansions of [Karney 2011](https://arxiv.org/abs/1002.1417), [Kawase 2011](http://www.gsi.go.jp/common/000062452.pdf), and [Kawase 2013](http://www.gsi.go.jp/common/000065826.pdf). The original DMATM implementaion are also present for reference.

## Bulding and testing

The code for the transformations is contained in the routines `utm.cpp` and `datum.cpp`. A test function is also included which can be compiled by seetting a flag.

To build the UTM test suite on Linux or MacOS use:

    g++ -DELLIPSE_TEST_MAIN -o ellipse_test utm.cpp datum.cpp -lm
    ./ellipse_test

This will output a series of foward and backward conversions between latitude/longitude and the UTM grid, reproducing table 2-11 of DMATM 8358.2 (which is included in this repository).
