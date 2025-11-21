#!/bin/bash

# Example usage script for Biosaur2 TOPP tool
# This script demonstrates how to use the Biosaur2 C++ implementation

set -e

echo "=========================================="
echo "Biosaur2 TOPP Tool - Example Usage"
echo "=========================================="
echo ""

# Check if executable exists
BIOSAUR2_EXE="./build/Biosaur2"

if [ ! -f "$BIOSAUR2_EXE" ]; then
    echo "ERROR: Biosaur2 executable not found!"
    echo "Please build the tool first using ./build.sh"
    exit 1
fi

echo "Biosaur2 executable found: $BIOSAUR2_EXE"
echo ""

# Show help
echo "=========================================="
echo "Displaying help information:"
echo "=========================================="
echo ""
$BIOSAUR2_EXE --help

echo ""
echo "=========================================="
echo "Example commands:"
echo "=========================================="
echo ""
echo "1. Basic usage:"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML"
echo ""
echo "2. With TSV output (compatible with Python biosaur2):"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -out_tsv output.features.tsv"
echo ""
echo "3. For short LC gradients (few minutes):"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -minlh 1"
echo ""
echo "4. For long LC gradients (60-180 min):"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -minlh 5"
echo ""
echo "5. For TOF data with higher mass tolerance:"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -htol 15.0 -itol 15.0"
echo ""
echo "6. Negative mode:"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -nm"
echo ""
echo "7. Custom m/z range and charge states:"
echo "   $BIOSAUR2_EXE -in input.mzML -out output.featureXML -minmz 400 -maxmz 1400 -cmin 2 -cmax 5"
echo ""

echo "=========================================="
echo "Parameter recommendations:"
echo "=========================================="
echo ""
echo "For 5-15 min LC gradients:"
echo "  -minlh 1-3"
echo ""
echo "For 60-180 min LC gradients:"
echo "  -minlh 5-10"
echo ""
echo "For high resolution data (Orbitrap):"
echo "  -htol 5.0 -itol 5.0"
echo ""
echo "For lower resolution data (TOF):"
echo "  -htol 10.0-15.0 -itol 10.0-15.0"
echo ""

echo "=========================================="
echo "To process a sample file:"
echo "=========================================="
echo ""
echo "Place your centroided mzML file in the current directory and run:"
echo "  $BIOSAUR2_EXE -in your_file.mzML -out your_file.featureXML -out_tsv your_file.features.tsv"
echo ""
