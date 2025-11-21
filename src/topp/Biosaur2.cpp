// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Mark Ivanov, Timo Sachsenberg $
// --------------------------------------------------------------------------

/**
 * @page TOPP_Biosaur2 Biosaur2
 *
 * @brief Feature detection for LC-MS1 data
 *
 * This TOPP tool is a C++ reimplementation of the Biosaur2 feature detection algorithm.
 * It detects peptide features in centroided LC-MS1 data by:
 * 1. Grouping peaks across scans into "hills"
 * 2. Splitting hills at valley points
 * 3. Detecting isotope patterns
 * 4. Calculating feature properties (m/z, RT, intensity, charge)
 *
 * Reference:
 * Abdrakhimov, et al. Biosaur: An open-source Python software for liquid 
 * chromatography-mass spectrometry peptide feature detection with ion mobility support.
 * https://doi.org/10.1002/rcm.9045
 *
 * <B>The command line parameters of this tool are:</B>
 * @verbinclude TOPP_Biosaur2.cli
 * <B>INI file documentation of this tool:</B>
 * @htmlinclude TOPP_Biosaur2.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/CoarseIsotopePatternGenerator.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <numeric>
#include <set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_Biosaur2 Biosaur2

  @brief Feature detection for LC-MS1 data using the Biosaur2 algorithm.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBiosaur2 :
  public TOPPBase
{
public:
  TOPPBiosaur2() :
    TOPPBase("Biosaur2", "Feature detection for LC-MS1 data", false)
  {
  }

protected:
  // Structure to represent a peak across multiple scans (a "hill")
  struct Hill
  {
    vector<Size> scan_indices;     // Indices of scans containing this hill
    vector<Size> peak_indices;     // Peak indices within each scan
    vector<double> mz_values;      // m/z values for each peak
    vector<double> intensities;    // Intensity values for each peak
    vector<double> rt_values;      // RT values for each scan
    double mz_median;              // Median m/z value
    double rt_start;               // Start RT
    double rt_end;                 // End RT
    double rt_apex;                // RT at apex
    double intensity_apex;         // Intensity at apex
    double intensity_sum;          // Sum of intensities
    Size length;                   // Number of scans
    Size hill_idx;                 // Unique hill identifier
  };

  // Structure to represent an isotope candidate
  struct IsotopeCandidate
  {
    Size hill_idx;                 // Hill index
    Size isotope_number;           // Isotope number (0=mono, 1=first 13C, etc)
    double mass_diff_ppm;          // Mass difference in ppm
    double cos_corr;               // Cosine correlation of RT profiles
  };

  // Structure to represent a feature
  struct PeptideFeature
  {
    double mz;                     // Monoisotopic m/z
    double rt_start;               // Start RT
    double rt_end;                 // End RT
    double rt_apex;                // RT at apex
    double intensity_apex;         // Intensity at apex
    double intensity_sum;          // Sum of intensities
    int charge;                    // Charge state
    Size n_isotopes;               // Number of isotopes
    Size n_scans;                  // Number of scans
    double mass_calib;             // Calibrated neutral mass
    vector<IsotopeCandidate> isotopes; // Isotope information
    Size mono_hill_idx;            // Monoisotopic hill index
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file (centroided data)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    
    registerOutputFile_("out", "<file>", "", "Output featureXML file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    
    registerOutputFile_("out_tsv", "<file>", "", "Optional: output TSV file (Biosaur2 format)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    // Parameters matching Python biosaur2
    registerDoubleOption_("mini", "<value>", 1.0, "Minimum intensity threshold", false);
    setMinFloat_("mini", 0.0);
    
    registerDoubleOption_("minmz", "<value>", 350.0, "Minimum m/z value", false);
    setMinFloat_("minmz", 0.0);
    
    registerDoubleOption_("maxmz", "<value>", 1500.0, "Maximum m/z value", false);
    setMinFloat_("maxmz", 0.0);
    
    registerDoubleOption_("htol", "<value>", 8.0, "Mass accuracy in ppm for combining peaks into hills", false);
    setMinFloat_("htol", 0.0);
    
    registerDoubleOption_("itol", "<value>", 8.0, "Mass accuracy in ppm for isotopic patterns", false);
    setMinFloat_("itol", 0.0);
    
    registerDoubleOption_("hvf", "<value>", 1.3, "Hill valley factor for splitting hills", false);
    setMinFloat_("hvf", 1.0);
    
    registerDoubleOption_("ivf", "<value>", 5.0, "Isotope valley factor for splitting isotope patterns (reserved for future use, not currently implemented)", false);
    setMinFloat_("ivf", 1.0);
    
    registerIntOption_("minlh", "<value>", 2, "Minimum number of scans for a hill", false);
    setMinInt_("minlh", 1);
    
    registerIntOption_("cmin", "<value>", 1, "Minimum charge state", false);
    setMinInt_("cmin", 1);
    
    registerIntOption_("cmax", "<value>", 6, "Maximum charge state", false);
    setMinInt_("cmax", 1);
    
    registerFlag_("nm", "Negative mode (affects neutral mass calculation)", false);
  }

  // Calculate ppm difference between two m/z values
  double calculatePPM(double mz1, double mz2)
  {
    return (mz1 - mz2) / mz2 * 1e6;
  }

  // Calculate median of a vector
  double calculateMedian(vector<double>& values)
  {
    if (values.empty()) return 0.0;
    
    vector<double> sorted = values;
    sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();
    
    if (n % 2 == 0)
    {
      return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    }
    else
    {
      return sorted[n/2];
    }
  }

  // Simple moving average filter
  vector<double> meanFilter(const vector<double>& data, Size window)
  {
    vector<double> result(data.size());
    Size half_window = window / 2;
    
    for (Size i = 0; i < data.size(); ++i)
    {
      Size start = (i >= half_window) ? i - half_window : 0;
      Size end = min(i + half_window + 1, data.size());
      
      double sum = 0.0;
      for (Size j = start; j < end; ++j)
      {
        sum += data[j];
      }
      result[i] = sum / (end - start);
    }
    
    return result;
  }

  // Cosine correlation between two RT profiles
  double cosineCorrelation(const vector<double>& intensities1, 
                          const vector<Size>& scans1,
                          const vector<double>& intensities2,
                          const vector<Size>& scans2)
  {
    // Build maps for efficient lookup
    map<Size, double> map1, map2;
    for (Size i = 0; i < scans1.size(); ++i)
    {
      map1[scans1[i]] = intensities1[i];
    }
    for (Size i = 0; i < scans2.size(); ++i)
    {
      map2[scans2[i]] = intensities2[i];
    }
    
    // Calculate cosine correlation
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;
    
    for (const auto& p1 : map1)
    {
      Size scan = p1.first;
      double i1 = p1.second;
      
      if (map2.find(scan) != map2.end())
      {
        double i2 = map2[scan];
        dot_product += i1 * i2;
      }
      norm1 += i1 * i1;
    }
    
    for (const auto& p2 : map2)
    {
      norm2 += p2.second * p2.second;
    }
    
    if (norm1 == 0.0 || norm2 == 0.0) return 0.0;
    
    return dot_product / (sqrt(norm1) * sqrt(norm2));
  }

  // Detect hills (peaks grouped across scans)
  vector<Hill> detectHills(const MSExperiment& exp, 
                          double htol_ppm,
                          double min_intensity,
                          double min_mz,
                          double max_mz)
  {
    vector<Hill> all_hills;
    map<Size, Hill> active_hills; // Map from hill_idx to hill
    Size next_hill_idx = 0;
    
    OPENMS_LOG_INFO << "Detecting hills across " << exp.size() << " MS1 spectra..." << endl;
    
    for (Size scan_idx = 0; scan_idx < exp.size(); ++scan_idx)
    {
      const MSSpectrum& spectrum = exp[scan_idx];
      
      // Only process MS1 spectra
      if (spectrum.getMSLevel() != 1) continue;
      
      double rt = spectrum.getRT();
      
      // Create new hills map for this scan
      map<Size, Hill> new_active_hills;
      set<Size> matched_peaks;
      
      // Try to extend existing hills
      for (auto& hill_pair : active_hills)
      {
        Hill& hill = hill_pair.second;
        double target_mz = hill.mz_median;
        double mz_tolerance = target_mz * htol_ppm * 1e-6;
        
        // Find closest peak within tolerance
        double best_diff = mz_tolerance + 1.0;
        Size best_peak_idx = 0;
        bool found = false;
        
        for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
        {
          // Skip peaks outside filter range
          double peak_mz = spectrum[peak_idx].getMZ();
          if (peak_mz < min_mz || peak_mz > max_mz) continue;
          if (spectrum[peak_idx].getIntensity() < min_intensity) continue;
          
          // Skip already matched peaks
          if (matched_peaks.find(peak_idx) != matched_peaks.end()) continue;
          
          double diff = abs(peak_mz - target_mz);
          if (diff <= mz_tolerance && diff < best_diff)
          {
            best_diff = diff;
            best_peak_idx = peak_idx;
            found = true;
          }
        }
        
        if (found)
        {
          // Extend hill
          matched_peaks.insert(best_peak_idx);
          hill.scan_indices.push_back(scan_idx);
          hill.peak_indices.push_back(best_peak_idx);
          hill.mz_values.push_back(spectrum[best_peak_idx].getMZ());
          hill.intensities.push_back(spectrum[best_peak_idx].getIntensity());
          hill.rt_values.push_back(rt);
          
          // Update median m/z
          vector<double> mz_copy = hill.mz_values;
          hill.mz_median = calculateMedian(mz_copy);
          
          new_active_hills[hill.hill_idx] = hill;
        }
        else
        {
          // Hill ended, save it
          all_hills.push_back(hill);
        }
      }
      
      // Start new hills from unmatched peaks
      for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
      {
        double peak_mz = spectrum[peak_idx].getMZ();
        double peak_int = spectrum[peak_idx].getIntensity();
        
        // Apply filters
        if (peak_mz < min_mz || peak_mz > max_mz) continue;
        if (peak_int < min_intensity) continue;
        if (matched_peaks.find(peak_idx) != matched_peaks.end()) continue;
        
        // Create new hill
        Hill new_hill;
        new_hill.hill_idx = next_hill_idx++;
        new_hill.scan_indices.push_back(scan_idx);
        new_hill.peak_indices.push_back(peak_idx);
        new_hill.mz_values.push_back(peak_mz);
        new_hill.intensities.push_back(peak_int);
        new_hill.rt_values.push_back(rt);
        new_hill.mz_median = peak_mz;
        
        new_active_hills[new_hill.hill_idx] = new_hill;
      }
      
      active_hills = new_active_hills;
    }
    
    // Save remaining active hills
    for (auto& hill_pair : active_hills)
    {
      all_hills.push_back(hill_pair.second);
    }
    
    OPENMS_LOG_INFO << "Detected " << all_hills.size() << " initial hills" << endl;
    
    return all_hills;
  }

  // Process hills: calculate properties and filter by minimum length
  vector<Hill> processHills(vector<Hill>& hills, Size min_length)
  {
    vector<Hill> processed_hills;
    
    for (auto& hill : hills)
    {
      hill.length = hill.scan_indices.size();
      
      // Filter by minimum length
      if (hill.length < min_length) continue;
      
      // Calculate RT properties
      hill.rt_start = hill.rt_values.front();
      hill.rt_end = hill.rt_values.back();
      
      // Find apex
      Size apex_idx = 0;
      double max_intensity = 0.0;
      for (Size i = 0; i < hill.intensities.size(); ++i)
      {
        if (hill.intensities[i] > max_intensity)
        {
          max_intensity = hill.intensities[i];
          apex_idx = i;
        }
      }
      hill.rt_apex = hill.rt_values[apex_idx];
      hill.intensity_apex = max_intensity;
      
      // Calculate sum intensity
      hill.intensity_sum = accumulate(hill.intensities.begin(), hill.intensities.end(), 0.0);
      
      processed_hills.push_back(hill);
    }
    
    OPENMS_LOG_INFO << "Processed " << processed_hills.size() << " hills (after filtering by min length)" << endl;
    
    return processed_hills;
  }

  // Split hills at valley points
  vector<Hill> splitHills(vector<Hill>& hills, double hvf, Size min_length)
  {
    vector<Hill> split_hills;
    Size next_hill_idx = 0;
    
    for (auto& hill : hills)
    {
      if (hill.length < min_length * 2)
      {
        // Too short to split, keep as is
        hill.hill_idx = next_hill_idx++;
        split_hills.push_back(hill);
        continue;
      }
      
      // Apply smoothing
      vector<double> smoothed = meanFilter(hill.intensities, 3);
      
      // Find local minima
      vector<Size> split_points;
      
      for (Size i = min_length - 1; i < hill.length - min_length; ++i)
      {
        // Check if this is a valley
        double valley_int = smoothed[i];
        if (valley_int == 0.0) continue;
        
        // Find max intensity on left and right
        double left_max = *max_element(smoothed.begin(), smoothed.begin() + i);
        double right_max = *max_element(smoothed.begin() + i + 1, smoothed.end());
        
        // Check valley criteria
        if (left_max / valley_int >= hvf && right_max / valley_int >= hvf)
        {
          // Check if not too close to previous split point
          if (split_points.empty() || i >= split_points.back() + min_length)
          {
            split_points.push_back(i);
          }
        }
      }
      
      if (split_points.empty())
      {
        // No split needed
        hill.hill_idx = next_hill_idx++;
        split_hills.push_back(hill);
      }
      else
      {
        // Split into multiple hills
        Size start_idx = 0;
        for (Size split_idx : split_points)
        {
          Hill new_hill;
          new_hill.hill_idx = next_hill_idx++;
          
          // Copy data from start to split point
          for (Size i = start_idx; i <= split_idx; ++i)
          {
            new_hill.scan_indices.push_back(hill.scan_indices[i]);
            new_hill.peak_indices.push_back(hill.peak_indices[i]);
            new_hill.mz_values.push_back(hill.mz_values[i]);
            new_hill.intensities.push_back(hill.intensities[i]);
            new_hill.rt_values.push_back(hill.rt_values[i]);
          }
          
          // Recalculate properties
          vector<double> mz_copy = new_hill.mz_values;
          new_hill.mz_median = calculateMedian(mz_copy);
          new_hill.length = new_hill.scan_indices.size();
          
          if (new_hill.length >= min_length)
          {
            split_hills.push_back(new_hill);
          }
          
          start_idx = split_idx + 1;
        }
        
        // Add remaining part
        if (start_idx < hill.length)
        {
          Hill new_hill;
          new_hill.hill_idx = next_hill_idx++;
          
          for (Size i = start_idx; i < hill.length; ++i)
          {
            new_hill.scan_indices.push_back(hill.scan_indices[i]);
            new_hill.peak_indices.push_back(hill.peak_indices[i]);
            new_hill.mz_values.push_back(hill.mz_values[i]);
            new_hill.intensities.push_back(hill.intensities[i]);
            new_hill.rt_values.push_back(hill.rt_values[i]);
          }
          
          vector<double> mz_copy = new_hill.mz_values;
          new_hill.mz_median = calculateMedian(mz_copy);
          new_hill.length = new_hill.scan_indices.size();
          
          if (new_hill.length >= min_length)
          {
            split_hills.push_back(new_hill);
          }
        }
      }
    }
    
    // Recalculate properties for all split hills
    vector<Hill> final_hills = processHills(split_hills, min_length);
    
    OPENMS_LOG_INFO << "After splitting: " << final_hills.size() << " hills" << endl;
    
    return final_hills;
  }

  // Detect isotope patterns and create features
  vector<PeptideFeature> detectIsotopePatterns(vector<Hill>& hills,
                                               double itol_ppm,
                                               int min_charge,
                                               int max_charge,
                                               bool negative_mode)
  {
    vector<PeptideFeature> features;
    set<Size> used_hills;
    
    OPENMS_LOG_INFO << "Detecting isotope patterns..." << endl;
    
    // Sort hills by m/z for efficient searching
    sort(hills.begin(), hills.end(), 
         [](const Hill& a, const Hill& b) { return a.mz_median < b.mz_median; });
    
    const double NEUTRON_MASS = 1.00286864;
    
    for (Size i = 0; i < hills.size(); ++i)
    {
      // Skip if already used
      if (used_hills.find(hills[i].hill_idx) != used_hills.end()) continue;
      
      const Hill& mono_hill = hills[i];
      double mono_mz = mono_hill.mz_median;
      
      // Try different charge states
      for (int charge = max_charge; charge >= min_charge; --charge)
      {
        double mz_spacing = NEUTRON_MASS / charge;
        
        vector<IsotopeCandidate> isotopes;
        bool pattern_valid = true;
        
        // Look for isotopes
        for (int iso_num = 1; iso_num <= 9; ++iso_num)
        {
          double expected_mz = mono_mz + iso_num * mz_spacing;
          double mz_tolerance = expected_mz * itol_ppm * 1e-6;
          
          // Search for matching hill
          bool found = false;
          for (Size j = i + 1; j < hills.size(); ++j)
          {
            // Stop if we're too far
            if (hills[j].mz_median > expected_mz + mz_tolerance) break;
            
            // Skip if already used
            if (used_hills.find(hills[j].hill_idx) != used_hills.end()) continue;
            
            double diff = abs(hills[j].mz_median - expected_mz);
            if (diff <= mz_tolerance)
            {
              // Check RT overlap
              double cos_corr = cosineCorrelation(
                mono_hill.intensities, mono_hill.scan_indices,
                hills[j].intensities, hills[j].scan_indices
              );
              
              if (cos_corr >= 0.6)
              {
                IsotopeCandidate candidate;
                candidate.hill_idx = hills[j].hill_idx;
                candidate.isotope_number = iso_num;
                candidate.mass_diff_ppm = calculatePPM(hills[j].mz_median, expected_mz);
                candidate.cos_corr = cos_corr;
                
                isotopes.push_back(candidate);
                found = true;
                break;
              }
            }
          }
          
          if (!found && iso_num == 1)
          {
            // Must have at least first isotope
            pattern_valid = false;
            break;
          }
          else if (!found)
          {
            // No more isotopes found
            break;
          }
        }
        
        // Create feature if we have a valid pattern
        if (pattern_valid && !isotopes.empty())
        {
          PeptideFeature feature;
          feature.mz = mono_mz;
          feature.rt_start = mono_hill.rt_start;
          feature.rt_end = mono_hill.rt_end;
          feature.rt_apex = mono_hill.rt_apex;
          feature.intensity_apex = mono_hill.intensity_apex;
          feature.intensity_sum = mono_hill.intensity_sum;
          feature.charge = charge;
          feature.n_isotopes = isotopes.size() + 1; // +1 for monoisotope
          feature.n_scans = mono_hill.length;
          feature.isotopes = isotopes;
          feature.mono_hill_idx = mono_hill.hill_idx;
          
          // Calculate neutral mass
          double proton_mass = 1.007276;
          if (negative_mode)
          {
            feature.mass_calib = mono_mz * charge + proton_mass * charge;
          }
          else
          {
            feature.mass_calib = mono_mz * charge - proton_mass * charge;
          }
          
          features.push_back(feature);
          
          // Mark hills as used
          used_hills.insert(mono_hill.hill_idx);
          for (const auto& iso : isotopes)
          {
            used_hills.insert(iso.hill_idx);
          }
          
          break; // Found valid charge state, move to next monoisotopic peak
        }
      }
    }
    
    OPENMS_LOG_INFO << "Detected " << features.size() << " features with isotope patterns" << endl;
    
    return features;
  }

  // Write features to TSV file (Biosaur2 format)
  void writeTSV(const vector<PeptideFeature>& features, const String& filename)
  {
    ofstream out(filename);
    
    // Write header
    out << "massCalib\trtApex\tintensityApex\tintensitySum\tcharge\t"
        << "nIsotopes\tnScans\tmz\trtStart\trtEnd" << endl;
    
    // Write features
    for (const auto& f : features)
    {
      out << f.mass_calib << "\t"
          << f.rt_apex << "\t"
          << f.intensity_apex << "\t"
          << f.intensity_sum << "\t"
          << f.charge << "\t"
          << f.n_isotopes << "\t"
          << f.n_scans << "\t"
          << f.mz << "\t"
          << f.rt_start << "\t"
          << f.rt_end << endl;
    }
    
    out.close();
    OPENMS_LOG_INFO << "Wrote " << features.size() << " features to TSV file: " << filename << endl;
  }

  // Convert features to OpenMS FeatureMap
  FeatureMap convertToFeatureMap(const vector<PeptideFeature>& features, const MSExperiment& exp)
  {
    FeatureMap feature_map;
    
    for (const auto& f : features)
    {
      Feature feature;
      feature.setMZ(f.mz);
      feature.setRT(f.rt_apex);
      feature.setIntensity(f.intensity_apex);
      feature.setCharge(f.charge);
      
      // Set quality (use number of isotopes as quality metric)
      feature.setOverallQuality(f.n_isotopes);
      
      // Create convex hull
      ConvexHull2D hull;
      vector<DPosition<2>> hull_points;
      hull_points.push_back(DPosition<2>(f.rt_start, f.mz));
      hull_points.push_back(DPosition<2>(f.rt_end, f.mz));
      hull.setHullPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
      
      // Add metadata
      feature.setMetaValue("mass_calib", f.mass_calib);
      feature.setMetaValue("n_isotopes", f.n_isotopes);
      feature.setMetaValue("n_scans", f.n_scans);
      feature.setMetaValue("intensity_sum", f.intensity_sum);
      
      feature_map.push_back(feature);
    }
    
    // Set protein identifications (empty)
    feature_map.getProteinIdentifications().resize(1);
    
    return feature_map;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_tsv = getStringOption_("out_tsv");
    
    double mini = getDoubleOption_("mini");
    double minmz = getDoubleOption_("minmz");
    double maxmz = getDoubleOption_("maxmz");
    double htol = getDoubleOption_("htol");
    double itol = getDoubleOption_("itol");
    double hvf = getDoubleOption_("hvf");
    // Note: ivf parameter is registered for API compatibility but not currently used
    // It would be used for splitting isotope patterns at local minima (future enhancement)
    // double ivf = getDoubleOption_("ivf");
    Size minlh = getIntOption_("minlh");
    int cmin = getIntOption_("cmin");
    int cmax = getIntOption_("cmax");
    bool negative_mode = getFlag_("nm");
    
    //-------------------------------------------------------------
    // Reading input
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading input file: " << in << endl;
    
    MSExperiment exp;
    MzMLFile mzml_file;
    mzml_file.load(in, exp);
    
    // Filter to MS1 only
    exp.getSpectra().erase(
      remove_if(exp.begin(), exp.end(), 
                [](const MSSpectrum& s) { return s.getMSLevel() != 1; }),
      exp.end()
    );
    
    OPENMS_LOG_INFO << "Loaded " << exp.size() << " MS1 spectra" << endl;
    
    if (exp.empty())
    {
      OPENMS_LOG_ERROR << "No MS1 spectra found in input file!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    
    //-------------------------------------------------------------
    // Feature detection
    //-------------------------------------------------------------
    
    // Step 1: Detect hills
    vector<Hill> hills = detectHills(exp, htol, mini, minmz, maxmz);
    
    // Step 2: Process hills (calculate properties, filter by length)
    hills = processHills(hills, minlh);
    
    // Step 3: Split hills at valleys
    hills = splitHills(hills, hvf, minlh);
    
    // Step 4: Detect isotope patterns
    vector<PeptideFeature> features = detectIsotopePatterns(hills, itol, cmin, cmax, negative_mode);
    
    //-------------------------------------------------------------
    // Writing output
    //-------------------------------------------------------------
    
    // Write FeatureXML
    FeatureMap feature_map = convertToFeatureMap(features, exp);
    FeatureXMLFile feature_file;
    feature_file.store(out, feature_map);
    OPENMS_LOG_INFO << "Wrote " << features.size() << " features to: " << out << endl;
    
    // Write TSV if requested
    if (!out_tsv.empty())
    {
      writeTSV(features, out_tsv);
    }
    
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPBiosaur2 tool;
  return tool.main(argc, argv);
}

/// @endcond
