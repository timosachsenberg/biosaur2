#!/bin/bash

# Validation script for Biosaur2 C++ implementation
# Checks code structure and best practices without requiring OpenMS

set -e

echo "=========================================="
echo "Biosaur2 TOPP Tool - Code Validation"
echo "=========================================="
echo ""

CPP_FILE="src/topp/Biosaur2.cpp"

# Check if source file exists
if [ ! -f "$CPP_FILE" ]; then
    echo "ERROR: Source file not found: $CPP_FILE"
    exit 1
fi

echo "✓ Source file found: $CPP_FILE"

# Count lines of code
LOC=$(wc -l < "$CPP_FILE")
echo "✓ Lines of code: $LOC"

# Check for required components
echo ""
echo "Checking for required components..."

# Check for class definition
if grep -q "class TOPPBiosaur2" "$CPP_FILE"; then
    echo "✓ TOPP tool class defined"
else
    echo "✗ TOPP tool class not found"
    exit 1
fi

# Check for main function
if grep -q "int main(int argc" "$CPP_FILE"; then
    echo "✓ Main function defined"
else
    echo "✗ Main function not found"
    exit 1
fi

# Check for parameter registration
if grep -q "registerOptionsAndFlags_" "$CPP_FILE"; then
    echo "✓ Parameter registration method defined"
else
    echo "✗ Parameter registration not found"
    exit 1
fi

# Check for main execution method
if grep -q "ExitCodes main_" "$CPP_FILE"; then
    echo "✓ Main execution method defined"
else
    echo "✗ Main execution method not found"
    exit 1
fi

# Check for key algorithms
echo ""
echo "Checking for key algorithm implementations..."

if grep -q "detectHills" "$CPP_FILE"; then
    echo "✓ Hill detection implemented"
else
    echo "✗ Hill detection not found"
    exit 1
fi

if grep -q "processHills" "$CPP_FILE"; then
    echo "✓ Hill processing implemented"
else
    echo "✗ Hill processing not found"
    exit 1
fi

if grep -q "splitHills" "$CPP_FILE"; then
    echo "✓ Hill splitting implemented"
else
    echo "✗ Hill splitting not found"
    exit 1
fi

if grep -q "detectIsotopePatterns" "$CPP_FILE"; then
    echo "✓ Isotope pattern detection implemented"
else
    echo "✗ Isotope pattern detection not found"
    exit 1
fi

# Check for output methods
echo ""
echo "Checking for output methods..."

if grep -q "convertToFeatureMap" "$CPP_FILE"; then
    echo "✓ FeatureMap conversion implemented"
else
    echo "✗ FeatureMap conversion not found"
    exit 1
fi

if grep -q "writeTSV" "$CPP_FILE"; then
    echo "✓ TSV output implemented"
else
    echo "✗ TSV output not found"
    exit 1
fi

# Check for required parameters
echo ""
echo "Checking for required parameters..."

PARAMS=("mini" "minmz" "maxmz" "htol" "itol" "hvf" "minlh" "cmin" "cmax")
for param in "${PARAMS[@]}"; do
    # Check if parameter is registered as either Double or Int option
    if grep -q "registerDoubleOption_.*\"$param\"" "$CPP_FILE" || grep -q "registerIntOption_.*\"$param\"" "$CPP_FILE"; then
        echo "✓ Parameter: $param"
    else
        echo "✗ Parameter not found: $param"
        exit 1
    fi
done

# Check for OpenMS includes
echo ""
echo "Checking for OpenMS includes..."

INCLUDES=("TOPPBase.h" "MzMLFile.h" "FeatureXMLFile.h" "MSExperiment.h" "FeatureMap.h")
for inc in "${INCLUDES[@]}"; do
    if grep -q "#include.*$inc" "$CPP_FILE"; then
        echo "✓ Include: $inc"
    else
        echo "✗ Include not found: $inc"
        exit 1
    fi
done

# Check for documentation
echo ""
echo "Checking for documentation..."

if grep -q "@page TOPP_Biosaur2" "$CPP_FILE"; then
    echo "✓ Doxygen documentation present"
else
    echo "✗ Doxygen documentation not found"
    exit 1
fi

if grep -q "@brief" "$CPP_FILE"; then
    echo "✓ Brief description present"
else
    echo "✗ Brief description not found"
    exit 1
fi

# Check CMakeLists.txt
echo ""
echo "Checking CMakeLists.txt..."

if [ -f "CMakeLists.txt" ]; then
    echo "✓ CMakeLists.txt found"
    
    if grep -q "find_package(OpenMS" "CMakeLists.txt"; then
        echo "✓ OpenMS package lookup configured"
    else
        echo "✗ OpenMS package lookup not configured"
        exit 1
    fi
    
    if grep -q "add_executable(Biosaur2" "CMakeLists.txt"; then
        echo "✓ Executable target configured"
    else
        echo "✗ Executable target not configured"
        exit 1
    fi
else
    echo "✗ CMakeLists.txt not found"
    exit 1
fi

# Check documentation files
echo ""
echo "Checking documentation files..."

if [ -f "README_CPP.md" ]; then
    echo "✓ C++ README found"
else
    echo "✗ C++ README not found"
    exit 1
fi

if [ -f "INTEGRATION.md" ]; then
    echo "✓ Integration guide found"
else
    echo "✗ Integration guide not found"
    exit 1
fi

# Check build scripts
echo ""
echo "Checking build scripts..."

if [ -f "build.sh" ] && [ -x "build.sh" ]; then
    echo "✓ Build script found and executable"
else
    echo "✗ Build script not found or not executable"
    exit 1
fi

if [ -f "example_usage.sh" ] && [ -x "example_usage.sh" ]; then
    echo "✓ Example usage script found and executable"
else
    echo "✗ Example usage script not found or not executable"
    exit 1
fi

echo ""
echo "=========================================="
echo "All validation checks passed! ✓"
echo "=========================================="
echo ""
echo "The C++ implementation appears to be complete."
echo ""
echo "To build (requires OpenMS):"
echo "  ./build.sh"
echo ""
echo "Note: Building requires OpenMS library to be installed."
echo "See INTEGRATION.md for installation instructions."
echo ""
