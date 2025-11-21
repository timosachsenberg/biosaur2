# Biosaur2 TOPP Tool - C++ Implementation

This directory contains a C++ reimplementation of Biosaur2 as an OpenMS TOPP tool.

## Overview

Biosaur2 is a feature detection algorithm for LC-MS1 data that:
1. Groups peaks across scans into "hills"
2. Splits hills at valley points
3. Detects isotope patterns
4. Calculates feature properties (m/z, RT, intensity, charge)

## Requirements

- OpenMS library (version 2.8 or higher)
- CMake (version 3.12 or higher)
- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)

## Building

### Using CMake

```bash
mkdir build
cd build
cmake ..
make
```

### Using OpenMS build system

If you want to integrate this into the OpenMS TOPP tools:

1. Copy `Biosaur2.cpp` to the OpenMS source tree under `src/topp/`
2. Rebuild OpenMS

```bash
cd /path/to/openms
mkdir build
cd build
cmake ..
make Biosaur2
```

## Usage

### Basic usage

```bash
Biosaur2 -in input.mzML -out output.featureXML
```

### With TSV output (compatible with Python biosaur2)

```bash
Biosaur2 -in input.mzML -out output.featureXML -out_tsv output.features.tsv
```

### Advanced usage

```bash
Biosaur2 -in input.mzML -out output.featureXML \
  -mini 10.0 \
  -minmz 400 \
  -maxmz 1400 \
  -htol 10.0 \
  -itol 10.0 \
  -minlh 3 \
  -cmin 2 \
  -cmax 5
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-in` | string | required | Input mzML file (centroided data) |
| `-out` | string | required | Output featureXML file |
| `-out_tsv` | string | optional | Output TSV file (Biosaur2 format) |
| `-mini` | float | 1.0 | Minimum intensity threshold |
| `-minmz` | float | 350.0 | Minimum m/z value |
| `-maxmz` | float | 1500.0 | Maximum m/z value |
| `-htol` | float | 8.0 | Mass accuracy in ppm for hills |
| `-itol` | float | 8.0 | Mass accuracy in ppm for isotopes |
| `-hvf` | float | 1.3 | Hill valley factor for splitting |
| `-ivf` | float | 5.0 | Isotope valley factor (reserved, not implemented) |
| `-minlh` | int | 2 | Minimum number of scans per hill |
| `-cmin` | int | 1 | Minimum charge state |
| `-cmax` | int | 6 | Maximum charge state |
| `-nm` | flag | false | Negative mode |

## Algorithm Overview

### 1. Hill Detection
- Peaks are grouped across consecutive scans based on m/z tolerance (`htol`)
- Each "hill" represents a series of peaks that track together across scans

### 2. Hill Processing
- Hills are filtered by minimum length (`minlh`)
- Apex intensity and RT are calculated
- Total intensity is summed

### 3. Hill Splitting
- Hills are split at valley points where the intensity drops significantly
- Valley detection uses the hill valley factor (`hvf`)
- Smoothing is applied before valley detection

### 4. Isotope Pattern Detection
- For each hill, the algorithm searches for isotope peaks
- Expected m/z spacing is calculated for each charge state
- RT profiles are compared using cosine correlation
- Features require at least the first isotope peak (M+1)

### 5. Feature Output
- Features are written in OpenMS FeatureXML format
- Optionally, TSV output compatible with Python biosaur2

## Differences from Python Implementation

This C++ implementation provides the core functionality of biosaur2 with some differences:

**Implemented:**
- Hill detection with mass tolerance
- Hill splitting at valleys
- Isotope pattern recognition
- Charge state determination
- Basic feature properties

**Not implemented (from Python version):**
- Ion mobility support (PASEF data)
- FAIMS support
- TOF-specific processing
- Profile mode processing
- Multiprocessing (C++ version uses single-threaded processing)
- Advanced isotope mass calibration

## Output Format

### FeatureXML (OpenMS format)
Standard OpenMS feature format with:
- m/z, RT, intensity, charge
- Convex hulls
- Metadata: mass_calib, n_isotopes, n_scans, intensity_sum

### TSV (Biosaur2 format)
Tab-separated file with columns:
- massCalib: Calibrated neutral mass
- rtApex: Retention time at apex
- intensityApex: Intensity at apex
- intensitySum: Sum of intensities
- charge: Charge state
- nIsotopes: Number of isotopes
- nScans: Number of scans
- mz: Monoisotopic m/z
- rtStart: Start retention time
- rtEnd: End retention time

## Performance

The C++ implementation is designed to be faster than the Python version:
- Native compiled code vs interpreted Python
- Efficient data structures
- No multiprocessing overhead

Expected performance:
- Small files (<1000 spectra): ~1-5 seconds
- Medium files (1000-10000 spectra): ~10-60 seconds
- Large files (>10000 spectra): ~1-10 minutes

## Reference

Abdrakhimov, et al. "Biosaur: An open-source Python software for liquid chromatography-mass spectrometry peptide feature detection with ion mobility support."
Rapid Communications in Mass Spectrometry, 2022.
https://doi.org/10.1002/rcm.9045

## Original Python Implementation

The original Python implementation is available at:
- https://github.com/markmipt/biosaur2

## License

This implementation follows the OpenMS license (3-clause BSD).
The original Biosaur2 is licensed under Apache 2.0.

## Author

- Implementation: Timo Sachsenberg
- Original Algorithm: Mark Ivanov

## Support

For issues with the C++ implementation, please open an issue in the repository.
For questions about the algorithm, refer to the original Biosaur2 publication.
