# Integrating Biosaur2 into OpenMS

This document describes how to integrate the Biosaur2 TOPP tool into the OpenMS framework.

## Option 1: Standalone Build (Recommended for Testing)

Build Biosaur2 as a standalone tool that links against an installed OpenMS:

```bash
# Make sure OpenMS is installed and available
# On Ubuntu/Debian:
sudo apt-get install libopenms-dev

# On macOS with Homebrew:
brew install openms

# Build Biosaur2
./build.sh

# The executable will be in build/Biosaur2
./build/Biosaur2 --help
```

## Option 2: Integration into OpenMS Source Tree

To include Biosaur2 as an official OpenMS TOPP tool:

### Step 1: Copy Source File

Copy `src/topp/Biosaur2.cpp` to the OpenMS source tree:

```bash
# Assuming you have OpenMS source at /path/to/openms
cp src/topp/Biosaur2.cpp /path/to/openms/src/topp/
```

### Step 2: Register in OpenMS Build

Edit `/path/to/openms/src/topp/CMakeLists.txt` and add Biosaur2 to the list of TOPP tools:

```cmake
# Find the section that lists TOPP executables
set(TOPP_executables
    # ... existing tools ...
    Biosaur2
    # ... more tools ...
)
```

### Step 3: Build OpenMS

```bash
cd /path/to/openms
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make Biosaur2
```

### Step 4: Install (Optional)

```bash
sudo make install
```

The Biosaur2 tool will now be available alongside other OpenMS TOPP tools.

## Option 3: Using OpenMS Development Environment

If you're developing within the OpenMS ecosystem:

### Using the OpenMS Docker Container

```bash
# Pull the OpenMS development image
docker pull openms/openms-dev:latest

# Run container with your code mounted
docker run -it -v /path/to/biosaur2:/workspace openms/openms-dev:latest

# Inside container, build Biosaur2
cd /workspace
mkdir build && cd build
cmake ..
make
```

### Using Conda Environment

```bash
# Create conda environment with OpenMS
conda create -n openms-dev -c bioconda openms

# Activate environment
conda activate openms-dev

# Build Biosaur2
./build.sh
```

## Verifying the Installation

After building, verify that Biosaur2 works correctly:

```bash
# Check if the tool runs
./build/Biosaur2 --help

# Expected output:
# Biosaur2 - Feature detection for LC-MS1 data
# Version: ...
# ...
```

## Usage in OpenMS Workflows

Once integrated, Biosaur2 can be used in OpenMS workflows:

### Standalone Usage

```bash
Biosaur2 -in data.mzML -out features.featureXML
```

### In a TOPP Workflow

```bash
# Example: Complete feature detection pipeline
FileFilter -in raw_data.mzML -out filtered.mzML -peak_options:level 1
Biosaur2 -in filtered.mzML -out features.featureXML
FeatureLinkerUnlabeledQT -in features.featureXML -out consensus.consensusXML
```

### Using KNIME

1. Install the OpenMS KNIME plugin
2. Drag the "Biosaur2" node into your workflow
3. Configure parameters in the node settings
4. Connect to upstream and downstream nodes

## Troubleshooting

### OpenMS Not Found

If CMake cannot find OpenMS:

```bash
# Set OpenMS_DIR environment variable
export OpenMS_DIR=/path/to/openms/share/OpenMS/cmake

# Or specify during cmake
cmake -DOpenMS_DIR=/path/to/openms/share/OpenMS/cmake ..
```

### Compilation Errors

Ensure you have:
- C++17 compatible compiler
- OpenMS version 2.8 or higher
- All OpenMS dependencies installed

### Runtime Errors

Check:
- Input file is centroided mzML
- Input file contains MS1 spectra
- OpenMS libraries are in your LD_LIBRARY_PATH (Linux) or DYLD_LIBRARY_PATH (macOS)

## Performance Tuning

For large datasets, consider:

```bash
# Use minimal parameters for faster processing
Biosaur2 -in large.mzML -out features.featureXML -minlh 3 -cmax 4

# Filter m/z range to reduce search space
Biosaur2 -in large.mzML -out features.featureXML -minmz 400 -maxmz 1200
```

## Contributing

If you want to contribute improvements to the Biosaur2 TOPP tool:

1. Make changes to `src/topp/Biosaur2.cpp`
2. Test thoroughly with various datasets
3. Update documentation in `README_CPP.md`
4. Submit a pull request

## Support

For issues with:
- **Biosaur2 C++ implementation**: Open issue in this repository
- **OpenMS integration**: OpenMS GitHub repository or mailing list
- **Algorithm questions**: See original Biosaur2 publication

## References

- OpenMS: https://www.openms.de/
- OpenMS Documentation: https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/index.html
- OpenMS GitHub: https://github.com/OpenMS/OpenMS
- Original Biosaur2: https://github.com/markmipt/biosaur2
