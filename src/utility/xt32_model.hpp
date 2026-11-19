/************************************************************************************************
  Sensor model of one Hesai PandarXT-32 unit (motiv patch), applied to every decoded point before it is published.

  Model xt32-v4 of convert_rings.py (the offline converter the calibration was measured with): this is a port of its
  correct_arrays(), operation by operation, so the driver and the converter give the same cloud to the float32 bit
  (checked on 742 M points). v4 has the parameters of v3 computed in two steps (11-19, user's choice: less work per
  point). Per point of ring r (0 = top):
    1. direction of the beam: firing rotation about z by the per-ring angle (cos / sin evaluated at load), and the
       global elevation offset as a rotation of (h, z) in the vertical plane of the point (cos / sin of the offset
       evaluated once at load: no atan2 / cos / sin per point);
    2. optical centre and near field together, in the frame of the point, u_h = (x, y) / h outwards and
       u_p = (-y, x) / h lateral:
       - optical centre (left, forward): the point is O + d t with |O + d t| = dist, i.e. R u_h + P u_p + k z e_z with
         R = k h + fwd, P = left. It keeps the measured range, so the near-field table is read at dist;
       - near field: p -= b p / |p| + a (-y, x, 0) / hypot(x, y), i.e. with f = 1 - b / dist, h' = sqrt(R^2 + P^2):
           A = f R + a P / h',  B = f P - a R / h',  z = f k z,  point = A u_h + B u_p + z e_z.
         Skipped where the table is zero to its end (dist >= r_zero), which is exact; f = 1 when b is all zero.
  x, y and z stay in double and are rounded to float32 once, at the end. Any change here must be made the same way in
  convert_rings.correct_arrays() (and checked bit for bit). Against convert_rings v3 (elevation applied last, about
  the origin) the points differ by 0.027 mm, almost all vertical (the optical centre sees the corrected elevation:
  31.5 mm * sin 0.05 deg).

  The file (driver/xt32_model_path) holds the calibration as in the scanner profile, block lidar_sensor_model
  (written by xt32_driver_model.py from convert_rings.py; motiv's lidar.launch.py copies it once to
  config/lidar/xt32_model_<serial>.yaml):
    model_version, lidar_hz, firing_scale, firing_da_rad ([] or 32 values), optical_centre_mm [left, forward],
    near_field_r_m, near_field_a_mm, near_field_b_mm, elev_offset_deg (and serial_number).
  The derived numbers are computed at load with the same operations as convert_rings.unit_params(), in the same
  order (firing angle -2 pi hz 1.512e-6 scale (r - 15.5) - da, mm * 1e-3, deg * (pi / 180), libm cos / sin), so
  they are the converter's doubles to the bit. The SDK must decode raw (no firetimes file, distance_correction_flag
  false, xt_spot_correction false, no transform), which the driver checks at start-up. Built without -ffast-math and
  without FMA (no -march), as the rest of the driver, so the arithmetic is the same IEEE double arithmetic numpy does.
************************************************************************************************/
#ifndef HESAI_ROS_XT32_MODEL_HPP_
#define HESAI_ROS_XT32_MODEL_HPP_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

class Xt32Model
{
public:
  // Loads the model; returns false and fills msg on any error (the model is then left disabled).
  bool Load(const std::string& path, std::string& msg)
  {
    enabled_ = false;
    double hz = 0.0, scale = 0.0, elev_deg = 0.0;
    std::vector<double> da, oc_mm, a_mm, b_mm;
    try {
      const YAML::Node m = YAML::LoadFile(path);
      version_ = m["model_version"].as<std::string>();
      serial_ = m["serial_number"].as<std::string>();
      hz = m["lidar_hz"].as<double>();
      scale = m["firing_scale"].as<double>();
      da = m["firing_da_rad"].as<std::vector<double>>();
      oc_mm = m["optical_centre_mm"].as<std::vector<double>>();
      r_ = m["near_field_r_m"].as<std::vector<double>>();
      a_mm = m["near_field_a_mm"].as<std::vector<double>>();
      b_mm = m["near_field_b_mm"].as<std::vector<double>>();
      elev_deg = m["elev_offset_deg"].as<double>();
    } catch (const std::exception& e) {
      msg = "xt32 model " + path + ": " + e.what();
      return false;
    }
    if (!da.empty() && da.size() != 32) {
      msg = "xt32 model " + path + ": firing_da_rad needs 0 or 32 values";
      return false;
    }
    if (oc_mm.size() != 2) {
      msg = "xt32 model " + path + ": optical_centre_mm needs [left, forward]";
      return false;
    }
    if (r_.size() != b_mm.size() || r_.size() != a_mm.size() || r_.size() == 1) {
      msg = "xt32 model " + path + ": near_field_r_m, _a_mm and _b_mm need the same length (0, or 2 or more)";
      return false;
    }
    // convert_rings.unit_params(), same operations in the same order
    const double base = -2.0 * M_PI * hz * 1.512e-6 * scale;
    cos_.assign(32, 1.0);
    sin_.assign(32, 0.0);
    for (int r = 0; r < 32; ++r) {
      double angle = base * (r - 15.5);
      if (!da.empty()) angle = angle - da[r];
      cos_[r] = std::cos(angle);
      sin_[r] = std::sin(angle);
    }
    left_ = oc_mm[0] * 1e-3;
    fwd_ = oc_mm[1] * 1e-3;
    a_.resize(a_mm.size());
    b_.resize(b_mm.size());
    for (size_t j = 0; j < r_.size(); ++j) {
      a_[j] = 1e-3 * a_mm[j];
      b_[j] = 1e-3 * b_mm[j];
    }
    elev_ = elev_deg * (M_PI / 180.0);
    std::ostringstream what;
    what << version_ << ", firing scale " << scale << (da.empty() ? "" : " + da") << ", optical centre " << oc_mm[0]
         << " / " << oc_mm[1] << " mm, near field " << r_.size() << " nodes, elevation offset " << elev_deg
         << " deg, " << hz << " Hz";
    model_ = what.str();
    for (size_t j = 1; j < r_.size(); ++j) {
      if (!(r_[j] > r_[j - 1])) {
        msg = "xt32 model " + path + ": near_field_r_m must be strictly ascending";
        return false;
      }
    }
    use_nf_ = !r_.empty();
    if (use_nf_) {
      // slopes as numpy.interp computes them; the table is zero from r_zero_ to its end (convert_rings.nf_zero_from)
      slope_a_.assign(r_.size(), 0.0);
      slope_b_.assign(r_.size(), 0.0);
      for (size_t j = 1; j < r_.size(); ++j) {
        slope_a_[j] = (a_[j] - a_[j - 1]) / (r_[j] - r_[j - 1]);
        slope_b_[j] = (b_[j] - b_[j - 1]) / (r_[j] - r_[j - 1]);
      }
      size_t j = r_.size();
      while (j > 0 && a_[j - 1] == 0.0 && b_[j - 1] == 0.0) --j;
      r_zero_ = j < r_.size() ? r_[j] : r_.back();
      has_b_ = std::any_of(b_.begin(), b_.end(), [](double v) { return v != 0.0; });
    }
    cos_e_ = std::cos(elev_);
    sin_e_ = std::sin(elev_);
    enabled_ = true;
    msg = "xt32 model of unit " + serial_ + " (" + path + "): " + model_ + "; two steps (firing + elevation, " +
          "optical centre + near field)";
    return true;
  }

  bool Enabled() const { return enabled_; }

  // Corrects one point in place (lidar frame, metres), as convert_rings.correct_arrays() (v4).
  void Apply(float& px, float& py, float& pz, uint16_t ring_in) const
  {
    const int ring = ring_in < 31 ? ring_in : 31;
    const double x = px, y = py;
    double z = pz;
    const double c = cos_[ring], s = sin_[ring];
    // 1. direction of the beam: firing rotation about z, then the elevation offset in the vertical plane
    double xr = c * x - s * y;
    double yr = s * x + c * y;
    double h = std::hypot(xr, yr);
    if (elev_ != 0.0 && h > 1e-6) {
      const double h2 = h * cos_e_ - z * sin_e_;
      const double z2 = h * sin_e_ + z * cos_e_;
      const double k = h2 / h;
      xr = xr * k;
      yr = yr * k;
      z = z2;
      h = h2;
    }
    // 2. optical centre + near field, in the frame of the point
    const double dist = std::sqrt(h * h + z * z);
    if (!(h > 1e-6)) {
      px = static_cast<float>(xr);
      py = static_cast<float>(yr);
      pz = static_cast<float>(z);
      return;
    }
    const double ux = xr / h, uy = yr / h;
    const bool oc = dist > 0.09 && std::isfinite(dist);
    double R = h, P = 0.0, zc = z, hp = h;
    if (oc) {
      const double o_t = fwd_ * h / dist;
      const double d = std::sqrt(o_t * o_t - (fwd_ * fwd_ + left_ * left_) + dist * dist) - o_t;
      const double k = d / dist;
      R = k * h + fwd_;
      P = left_;
      zc = k * z;
      hp = std::sqrt(R * R + P * P);
    }
    double A = R, B = P;
    bool moved = oc;
    if (use_nf_ && hp > 1e-6 && std::isfinite(dist) && dist < r_zero_) {
      double a, b = 0.0;
      Interp(dist, a, b);
      if (has_b_) {
        const double f = 1.0 - b / dist;
        A = f * R + a * P / hp;
        B = f * P - a * R / hp;
        zc = f * zc;
      } else {  // f = 1 exactly
        A = R + a * P / hp;
        B = P - a * R / hp;
      }
      moved = true;
    }
    if (moved) {
      xr = A * ux - B * uy;
      yr = A * uy + B * ux;
    }
    px = static_cast<float>(xr);
    py = static_cast<float>(yr);
    pz = static_cast<float>(zc);
  }

private:
  // numpy.interp for rho < last node: first value at or below the first node, exact node values, linear between.
  void Interp(double rho, double& a, double& b) const
  {
    if (rho <= r_.front()) {
      a = a_.front();
      b = b_.front();
      return;
    }
    // first node above rho among r_[1 .. n-2], else n-1 (same segment as the linear search of numpy.interp)
    const size_t j = std::upper_bound(r_.begin() + 1, r_.end() - 1, rho) - r_.begin();
    if (r_[j - 1] == rho) {
      a = a_[j - 1];
      b = b_[j - 1];
      return;
    }
    a = slope_a_[j] * (rho - r_[j - 1]) + a_[j - 1];
    b = has_b_ ? slope_b_[j] * (rho - r_[j - 1]) + b_[j - 1] : 0.0;
  }

  bool enabled_ = false, use_nf_ = false, has_b_ = false;
  std::vector<double> cos_, sin_, r_, b_, a_, slope_a_, slope_b_;
  double left_ = 0.0, fwd_ = 0.0, elev_ = 0.0, cos_e_ = 1.0, sin_e_ = 0.0, r_zero_ = 0.0;
  std::string model_, serial_, version_;
};

#endif  // HESAI_ROS_XT32_MODEL_HPP_
