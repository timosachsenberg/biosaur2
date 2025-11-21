# Biosaur2: Python vs C++ Implementation Comparison

This document compares the Python and C++ implementations of Biosaur2.

## Architecture Comparison

### Python Implementation
- **Language**: Python 3 with Cython optimizations
- **Dependencies**: NumPy, SciPy, pandas, pyteomics
- **Performance**: Uses multiprocessing for parallelization
- **File I/O**: pyteomics for mzML reading
- **Output**: TSV files

### C++ Implementation
- **Language**: C++17
- **Dependencies**: OpenMS library
- **Performance**: Native compiled code with OpenMP multithreading support
- **File I/O**: OpenMS MzMLFile
- **Output**: FeatureXML (OpenMS) and TSV formats

## Parameter Mapping

| Python Parameter | C++ Parameter | Description |
|-----------------|---------------|-------------|
| `-mini` | `-mini` | Minimum intensity threshold |
| `-minmz` | `-minmz` | Minimum m/z value |
| `-maxmz` | `-maxmz` | Maximum m/z value |
| `-htol` | `-htol` | Mass accuracy for hills (ppm) |
| `-itol` | `-itol` | Mass accuracy for isotopes (ppm) |
| `-hvf` | `-hvf` | Hill valley factor |
| `-ivf` | `-ivf` | Isotope valley factor (both implemented) |
| `-minlh` | `-minlh` | Minimum hill length |
| `-cmin` | `-cmin` | Minimum charge state |
| `-cmax` | `-cmax` | Maximum charge state |
| `-nm` | `-nm` | Negative mode |
| `-nprocs` | `-threads` | Number of threads/processes |
| `-tof` | `-tof` | TOF-specific processing |
| `-write_hills` | `-write_hills` | Hills output |
| `-iuse` | `-iuse` | Isotopes for intensity calc |
| `-ignore_iso_calib` | `-ignore_iso_calib` | Turn off isotope calibration |
| `-use_hill_calib` | `-use_hill_calib` | Turn on hill calibration |
| `-profile` | `-profile` | Profile mode processing |
| `-write_extra_details` | N/A | Extra feature details (Python only) |
| `-o` | `-out_tsv` | TSV output path |
| N/A | `-out` | FeatureXML output path |
| N/A | `-out_hills` | Hills TSV output path |

## Feature Support Matrix

| Feature | Python | C++ | Notes |
|---------|--------|-----|-------|
| **Core Algorithm** |
| Hill detection | ✓ | ✓ | Both implementations |
| Hill splitting | ✓ | ✓ | Valley-based splitting |
| Isotope detection | ✓ | ✓ | Charge state analysis |
| Feature calculation | ✓ | ✓ | RT, m/z, intensity |
| **Advanced Features** |
| Ion mobility (PASEF) | ✓ | ✓ | Both implementations (C++ uses FloatDataArray) |
| FAIMS support | ✓ | ✓ | Both implementations (C++ tracks drift time) |
| TOF processing | ✓ | ✓ | Both implementations |
| Profile mode | ✓ | ✓ | Both implementations (C++ uses PeakPickerHiRes) |
| Multiprocessing | ✓ | ✓ | Both (C++ uses OpenMP) |
| Isotope splitting (ivf) | ✓ | ✓ | Both implementations |
| Isotopes for intensity (iuse) | ✓ | ✓ | Both implementations |
| **Output Formats** |
| TSV | ✓ | ✓ | Both implementations |
| FeatureXML | ✗ | ✓ | C++ only |
| Hills output | ✓ | ✓ | Both implementations |
| **Calibration** |
| Hill mass calibration | ✓ | ✓ | Both implementations |
| Isotope mass calibration | ✓ | ✓ | Both implementations (C++ auto-enabled by default) |
| **Constants** |
| OpenMS constants | N/A | ✓ | C++ uses Constants::C13C12_MASSDIFF_U |

## Algorithm Mapping

### 1. Data Loading

**Python:**
```python
def process_mzml(args):
    for z in MS1OnlyMzML(source=input_mzml_path):
        # Filter by intensity, m/z range
        # Store in list
```

**C++:**
```cpp
MSExperiment exp;
MzMLFile mzml_file;
mzml_file.load(in, exp);
// Filter to MS1 only
```

### 2. Hill Detection

**Python (cutils.pyx):**
```python
def detect_hills(data_for_analyse_tmp, args, mz_step, paseftol):
    # Group peaks across scans using mz_step tolerance
    # Build hills dictionary with scan/peak tracking
```

**C++:**
```cpp
vector<Hill> detectHills(const MSExperiment& exp, 
                        double htol_ppm, ...)
{
    // Track active hills
    // Extend or close hills each scan
    // Return completed hills
}
```

### 3. Hill Processing

**Python:**
```python
def process_hills(hills_dict, data_for_analyse_tmp, 
                  mz_step, paseftol, args):
    # Calculate median m/z
    # Find apex intensity and RT
    # Calculate properties
```

**C++:**
```cpp
vector<Hill> processHills(vector<Hill>& hills, Size min_length)
{
    // Calculate RT properties
    // Find apex
    // Sum intensities
    // Filter by length
}
```

### 4. Hill Splitting

**Python (cutils.pyx):**
```python
def split_peaks(hills_dict, data_for_analyse_tmp, args, ...):
    # Apply smoothing filter
    # Find local minima (valleys)
    # Split at valley points
    # Check hvf threshold
```

**C++:**
```cpp
vector<Hill> splitHills(vector<Hill>& hills, double hvf, Size min_length)
{
    // Mean filter for smoothing
    // Detect valleys
    // Split into new hills
    // Recalculate properties
}
```

### 5. Isotope Pattern Detection

**Python (cutils.pyx):**
```python
def get_initial_isotopes(hills_dict, isotopes_mass_accuracy, ...):
    # For each hill (monoisotope candidate)
    # Search for isotopes at expected m/z
    # Calculate cosine correlation
    # Build isotope list
```

**C++:**
```cpp
vector<PeptideFeature> detectIsotopePatterns(vector<Hill>& hills, ...)
{
    // Try different charge states
    // Search for isotope peaks
    // Check RT correlation
    // Create features
}
```

## Data Structure Mapping

### Python Hills Dictionary
```python
hills_dict = {
    'hills_idx_array': [...],        # Hill indices
    'scan_idx_array': [...],         # Scan indices
    'orig_idx_array': [...],         # Peak indices
    'mzs_array': [...],              # m/z values
    'intensity_array': [...],        # Intensities
    'hills_mz_median': [...],        # Median m/z per hill
    # ... more fields
}
```

### C++ Hill Structure
```cpp
struct Hill {
    vector<Size> scan_indices;
    vector<Size> peak_indices;
    vector<double> mz_values;
    vector<double> intensities;
    vector<double> rt_values;
    double mz_median;
    double rt_start, rt_end, rt_apex;
    double intensity_apex, intensity_sum;
    Size length;
    Size hill_idx;
};
```

## Performance Comparison

### Python
- **Advantages**: 
  - Multiprocessing for large files
  - Optimized Cython code for critical sections
- **Speed**: Fast for small-medium files, scales well with CPU cores
- **Memory**: Can be high due to Python overhead

### C++
- **Advantages**:
  - Native compiled code
  - No Python/GIL overhead
  - Efficient memory usage
  - OpenMP multithreading support
- **Speed**: Very fast for all file sizes, scales with thread count
- **Memory**: Lower overhead, more efficient

### Estimated Performance

| Dataset Size | Python (4 cores) | C++ (1 core) | C++ (4 cores) |
|-------------|------------------|--------------|---------------|
| 500 spectra | ~5 sec | ~1-2 sec | ~0.5-1 sec |
| 5000 spectra | ~30 sec | ~10-20 sec | ~5-10 sec |
| 50000 spectra | ~5 min | ~2-5 min | ~1-2.5 min |

*Note: Actual performance depends on data complexity, m/z range, and parameters*

## Output Comparison

### Python TSV Output
```tsv
massCalib  rtApex  intensityApex  intensitySum  charge  nIsotopes  nScans  mz  rtStart  rtEnd  FAIMS  im  ...
```

### C++ TSV Output
```tsv
massCalib  rtApex  intensityApex  intensitySum  charge  nIsotopes  nScans  mz  rtStart  rtEnd
```

The C++ TSV output is compatible with Python output (subset of columns).

### C++ FeatureXML Output
OpenMS standard format:
- Compatible with downstream OpenMS tools
- Contains convex hulls
- Includes metadata
- Can be visualized in TOPPView

## Use Case Recommendations

### Use Python Implementation When:
- You need ion mobility support (PASEF data)
- You need FAIMS support
- You want to use the advanced calibration features
- You're working in a Python ecosystem
- You need hills output for analysis

### Use C++ Implementation When:
- You're using OpenMS workflows
- You need maximum performance
- You want OpenMS FeatureXML output
- You're building a C++ application
- You want native KNIME integration

## Migration Guide

### From Python to C++

Python command:
```bash
biosaur2 input.mzML -mini 10 -minlh 3 -htol 10 -o output.features.tsv
```

Equivalent C++ command:
```bash
Biosaur2 -in input.mzML -mini 10 -minlh 3 -htol 10 \
  -out output.featureXML -out_tsv output.features.tsv
```

### From C++ to Python

C++ command:
```bash
Biosaur2 -in input.mzML -out features.featureXML -out_tsv features.tsv
```

Equivalent Python command:
```bash
biosaur2 input.mzML -o features.tsv
```

*Note: To read FeatureXML in Python, you'll need to convert it first*

## Conclusion

Both implementations provide robust feature detection:
- **Python**: Full-featured with advanced options
- **C++**: High-performance core algorithm with OpenMS integration

Choose based on your workflow requirements and ecosystem.
