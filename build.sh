#!/bin/bash

# Build script for Biosaur2 TOPP tool
# This script helps build the Biosaur2 C++ TOPP tool

set -e

echo "=========================================="
echo "Biosaur2 TOPP Tool - Build Script"
echo "=========================================="
echo ""

# Check if OpenMS is available
if ! command -v openms-config &> /dev/null; then
    echo "Warning: openms-config not found in PATH"
    echo "Make sure OpenMS is installed and available"
    echo ""
fi

# Create build directory
BUILD_DIR="build"
if [ -d "$BUILD_DIR" ]; then
    echo "Removing existing build directory..."
    rm -rf "$BUILD_DIR"
fi

echo "Creating build directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Run CMake
echo ""
echo "Running CMake..."
cmake .. || {
    echo ""
    echo "ERROR: CMake configuration failed!"
    echo "Please ensure OpenMS is installed and the OpenMS_DIR is set correctly."
    echo ""
    echo "You can try setting OpenMS_DIR manually:"
    echo "  export OpenMS_DIR=/path/to/openms/share/OpenMS/cmake"
    echo "  or"
    echo "  cmake -DOpenMS_DIR=/path/to/openms/share/OpenMS/cmake .."
    echo ""
    exit 1
}

# Build
echo ""
echo "Building Biosaur2..."
make -j$(nproc 2>/dev/null || echo 2) || {
    echo ""
    echo "ERROR: Build failed!"
    exit 1
}

echo ""
echo "=========================================="
echo "Build completed successfully!"
echo "=========================================="
echo ""
echo "The Biosaur2 executable is located at:"
echo "  $BUILD_DIR/Biosaur2"
echo ""
echo "To install system-wide, run:"
echo "  sudo make install"
echo ""
echo "To test the tool, run:"
echo "  ./Biosaur2 --help"
echo ""
