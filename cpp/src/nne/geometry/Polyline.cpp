#include "nne/geometry/Polyline.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace fastnsa::nne {
namespace {

template <typename T>
T read_scalar(const std::uint8_t*& cursor, const std::uint8_t* end, bool little_endian) {
    static_assert(std::is_trivially_copyable_v<T>);
    if (static_cast<std::size_t>(end - cursor) < sizeof(T)) {
        throw std::invalid_argument("Truncated LineString WKB.");
    }
    std::uint8_t bytes[sizeof(T)];
    std::memcpy(bytes, cursor, sizeof(T));
    cursor += sizeof(T);

    const std::uint16_t marker = 1;
    const bool host_little_endian = *reinterpret_cast<const std::uint8_t*>(&marker) == 1;
    if (host_little_endian != little_endian) {
        std::reverse(bytes, bytes + sizeof(T));
    }

    T value;
    std::memcpy(&value, bytes, sizeof(T));
    return value;
}

double segment_length(const Point2D& first, const Point2D& second) noexcept {
    return std::hypot(second.x - first.x, second.y - first.y);
}

}  // namespace

double Envelope::distance_squared(const Point2D& point) const noexcept {
    const double dx = point.x < min_x ? min_x - point.x
        : point.x > max_x ? point.x - max_x : 0.0;
    const double dy = point.y < min_y ? min_y - point.y
        : point.y > max_y ? point.y - max_y : 0.0;
    return dx * dx + dy * dy;
}

void Envelope::expand(const Envelope& other) noexcept {
    min_x = std::min(min_x, other.min_x);
    min_y = std::min(min_y, other.min_y);
    max_x = std::max(max_x, other.max_x);
    max_y = std::max(max_y, other.max_y);
}

Polyline::Polyline(std::vector<Point2D> points) : points_(std::move(points)) {
    if (points_.size() < 2) {
        throw std::invalid_argument("A network edge must contain at least two vertices.");
    }
    for (const auto& point : points_) {
        if (!std::isfinite(point.x) || !std::isfinite(point.y)) {
            throw std::invalid_argument("Network edge coordinates must be finite.");
        }
    }

    envelope_ = {points_[0].x, points_[0].y, points_[0].x, points_[0].y};
    cumulative_lengths_.resize(points_.size(), 0.0);
    for (std::size_t index = 1; index < points_.size(); ++index) {
        envelope_.min_x = std::min(envelope_.min_x, points_[index].x);
        envelope_.min_y = std::min(envelope_.min_y, points_[index].y);
        envelope_.max_x = std::max(envelope_.max_x, points_[index].x);
        envelope_.max_y = std::max(envelope_.max_y, points_[index].y);
        length_ += segment_length(points_[index - 1], points_[index]);
        cumulative_lengths_[index] = length_;
    }
}

Projection Polyline::project(const Point2D& point) const noexcept {
    Projection best{std::numeric_limits<double>::infinity(), 0.0};
    for (std::size_t index = 1; index < points_.size(); ++index) {
        const auto& first = points_[index - 1];
        const auto& second = points_[index];
        const double dx = second.x - first.x;
        const double dy = second.y - first.y;
        const double length_squared = dx * dx + dy * dy;
        double fraction = 0.0;
        if (length_squared > 0.0) {
            fraction = ((point.x - first.x) * dx + (point.y - first.y) * dy) /
                length_squared;
            fraction = std::clamp(fraction, 0.0, 1.0);
        }
        const double projected_x = first.x + fraction * dx;
        const double projected_y = first.y + fraction * dy;
        const double distance_x = point.x - projected_x;
        const double distance_y = point.y - projected_y;
        const double distance_squared = distance_x * distance_x + distance_y * distance_y;
        if (distance_squared < best.distance_squared) {
            best.distance_squared = distance_squared;
            best.offset = cumulative_lengths_[index - 1] +
                fraction * (cumulative_lengths_[index] - cumulative_lengths_[index - 1]);
        }
    }
    best.offset = std::clamp(best.offset, 0.0, length_);
    return best;
}

Polyline parse_linestring_wkb(const std::uint8_t* data, std::size_t size) {
    if (!data || size < 9) {
        throw std::invalid_argument("LineString WKB is empty or truncated.");
    }
    const auto* cursor = data;
    const auto* end = data + size;
    const std::uint8_t byte_order = *cursor++;
    if (byte_order > 1) {
        throw std::invalid_argument("LineString WKB has an invalid byte order.");
    }
    const bool little_endian = byte_order == 1;
    const std::uint32_t encoded_type = read_scalar<std::uint32_t>(cursor, end, little_endian);

    const bool ewkb_z = (encoded_type & 0x80000000u) != 0;
    const bool ewkb_m = (encoded_type & 0x40000000u) != 0;
    const bool ewkb_srid = (encoded_type & 0x20000000u) != 0;
    std::uint32_t geometry_type = encoded_type & 0x1fffffffu;
    int dimensions = 2 + static_cast<int>(ewkb_z) + static_cast<int>(ewkb_m);
    if (geometry_type >= 3000) {
        geometry_type -= 3000;
        dimensions = 4;
    } else if (geometry_type >= 2000) {
        geometry_type -= 2000;
        dimensions = 3;
    } else if (geometry_type >= 1000) {
        geometry_type -= 1000;
        dimensions = 3;
    }
    if (geometry_type != 2) {
        throw std::invalid_argument("Only LineString WKB is supported for network edges.");
    }
    if (ewkb_srid) {
        (void)read_scalar<std::uint32_t>(cursor, end, little_endian);
    }

    const std::uint32_t point_count = read_scalar<std::uint32_t>(cursor, end, little_endian);
    if (point_count < 2) {
        throw std::invalid_argument("A LineString WKB edge must contain at least two vertices.");
    }
    const std::size_t required = static_cast<std::size_t>(point_count) *
        static_cast<std::size_t>(dimensions) * sizeof(double);
    if (static_cast<std::size_t>(end - cursor) < required) {
        throw std::invalid_argument("LineString WKB coordinate data is truncated.");
    }

    std::vector<Point2D> points;
    points.reserve(point_count);
    for (std::uint32_t index = 0; index < point_count; ++index) {
        const double x = read_scalar<double>(cursor, end, little_endian);
        const double y = read_scalar<double>(cursor, end, little_endian);
        for (int dimension = 2; dimension < dimensions; ++dimension) {
            (void)read_scalar<double>(cursor, end, little_endian);
        }
        points.push_back({x, y});
    }
    return Polyline(std::move(points));
}

}  // namespace fastnsa::nne
