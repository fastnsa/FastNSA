#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fastnsa::nne {

struct Point2D {
    double x = 0.0;
    double y = 0.0;
};

struct Envelope {
    double min_x = 0.0;
    double min_y = 0.0;
    double max_x = 0.0;
    double max_y = 0.0;

    double distance_squared(const Point2D& point) const noexcept;
    void expand(const Envelope& other) noexcept;
};

struct Projection {
    double distance_squared = 0.0;
    double offset = 0.0;
};

class Polyline {
public:
    explicit Polyline(std::vector<Point2D> points);

    const Envelope& envelope() const noexcept { return envelope_; }
    double length() const noexcept { return length_; }
    Projection project(const Point2D& point) const noexcept;

private:
    std::vector<Point2D> points_;
    std::vector<double> cumulative_lengths_;
    Envelope envelope_;
    double length_ = 0.0;
};

Polyline parse_linestring_wkb(const std::uint8_t* data, std::size_t size);

}  // namespace fastnsa::nne
