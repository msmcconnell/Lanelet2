#pragma once
#include <boost/geometry/algorithms/intersection.hpp>
#include "../../primitives/LineString.h"
#include "../../primitives/Traits.h"
#include "../GeometryHelper.h"
#include "../Point.h"

#include <iostream>

namespace lanelet {
namespace geometry {
namespace internal {
inline auto crossProd(const BasicPoint3d& p1, const BasicPoint3d& p2) { return p1.cross(p2).eval(); }
inline auto crossProd(const BasicPoint2d& p1, const BasicPoint2d& p2) {
  return BasicPoint3d(p1.x(), p1.y(), 0.).cross(BasicPoint3d(p2.x(), p2.y(), 0.)).eval();
}

template <typename LineStringT, typename BasicPointT>
auto findPoint(const LineStringT& ls, const BasicPointT& p) {
  return std::find_if(ls.begin(), ls.end(), [&p](const auto& elem) { return boost::geometry::equals(elem, p); });
}

template <typename PointT>
bool pointIsLeftOf(const PointT& pSeg1, const PointT& pSeg2, const PointT& p) {
  return crossProd(PointT(pSeg2 - pSeg1), PointT(p - pSeg1)).z() > 0;
}

template <typename LineStringT>
LineStringT invert(const LineStringT& ls) {
  return ls.invert();
}

template <>
inline BasicLineString2d invert(const BasicLineString2d& ls) {
  return BasicLineString2d{ls.rbegin(), ls.rend()};
}

template <>
inline BasicLineString3d invert(const BasicLineString3d& ls) {
  return BasicLineString3d{ls.rbegin(), ls.rend()};
}

template <typename LineStringT, typename BasicPointT>
bool isLeftOf(const LineStringT& ls, const BasicPointT& p, const helper::ProjectedPoint<BasicPointT>& projectedPoint) {
  BasicPointT pSeg1 = projectedPoint.result->segmentPoint1;
  BasicPointT pSeg2 = projectedPoint.result->segmentPoint2;
  BasicPointT projPoint = projectedPoint.result->projectedPoint;
  bool isLeft = pointIsLeftOf(pSeg1, pSeg2, p);
  if (pSeg2 == projPoint) {
    auto nextSegPointIt = std::next(findPoint(ls, pSeg2));
    if (nextSegPointIt != ls.end()) {
      // see stackoverflow.com/questions/10583212
      BasicPointT nextSegPoint;
      boost::geometry::convert(*nextSegPointIt, nextSegPoint);
      if (isLeft != pointIsLeftOf(pSeg2, nextSegPoint, p) && isLeft == pointIsLeftOf(pSeg1, pSeg2, nextSegPoint)) {
        return !isLeft;
      }
    }
  }
  return isLeft;
}

template <typename LineStringT, typename PointT>
std::pair<double, helper::ProjectedPoint<PointT>> signedDistanceImpl(const LineStringT lineString, const PointT& p) {
  using BasicPoint = PointT;
  helper::ProjectedPoint<BasicPoint> projectedPoint;
  const auto d = distance(lineString, p, projectedPoint);
  auto isLeft = isLeftOf(lineString, p, projectedPoint);
  return {isLeft ? d : -d, projectedPoint};
}

std::pair<BasicPoint3d, BasicPoint3d> projectedPoint3d(const ConstHybridLineString3d& l1,
                                                       const ConstHybridLineString3d& l2);

std::pair<BasicPoint3d, BasicPoint3d> projectedPoint3d(const CompoundHybridLineString3d& l1,
                                                       const CompoundHybridLineString3d& l2);

template <typename HybridLineStringT>
BasicPoint2d fromArcCoords(const HybridLineStringT& hLineString, const BasicPoint2d& projStart, const int startIdx,
                           const int endIdx, const double distance) {
  const auto dx(hLineString[endIdx](0) - hLineString[startIdx](0));
  const auto dy(hLineString[endIdx](1) - hLineString[startIdx](1));

  return projStart + Eigen::Vector2d(-dy, dx).normalized() * distance;
}

template <typename LineString2dT>
BasicPoint2d shiftLateral(const LineString2dT& lineString, const int idx, const double offset) {
  Eigen::Vector2d perpendicular;
  double realOffset = offset;
  if (idx == 0)
    perpendicular = utils::toBasicPoint(lineString[idx + 1]) - utils::toBasicPoint(lineString[idx]);
  else if (idx == lineString.size() - 1)
    perpendicular = utils::toBasicPoint(lineString[idx]) - utils::toBasicPoint(lineString[idx - 1]);
  else {
    Eigen::Vector2d following =
        (utils::toBasicPoint(lineString[idx + 1]) - utils::toBasicPoint(lineString[idx])).normalized();
    Eigen::Vector2d preceding =
        (utils::toBasicPoint(lineString[idx]) - utils::toBasicPoint(lineString[idx - 1])).normalized();
    perpendicular = following + preceding;
    realOffset = offset / std::cos(0.5 * std::acos(following.dot(preceding)));
  }

  Eigen::Vector2d direction(-perpendicular(1), perpendicular(0));
  return utils::toBasicPoint(lineString[idx]) + direction.normalized() * realOffset;
}

template <typename LineStringT>
std::vector<size_t> sortAlongSIdxs(const LineStringT& ls, const BasicPoints2d& points) {
  if (points.empty()) return std::vector<size_t>();
  std::vector<std::pair<size_t, double>> projections(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    const auto& p = points.at(i);
    auto arcC = toArcCoordinates(ls, p);
    projections[i] = std::make_pair(i, arcC.length);
  }
  std::sort(projections.begin(), projections.end(),
            [](const std::pair<size_t, double>& lhs, const std::pair<size_t, double>& rhs) -> bool {
              return lhs.second < rhs.second;
            });
  std::vector<size_t> ret;
  for (const auto& p : projections) {
    ret.emplace_back(p.first);
  }
  return ret;
}

template <typename LineStringT>
BasicPoints2d sortAlongS(const LineStringT& ls, const BasicPoints2d& points) {
  auto idxs = sortAlongSIdxs(ls, points);
  BasicPoints2d ret;
  for (const auto& i : idxs) ret.emplace_back(points.at(i));
  return ret;
}

}  // namespace internal

template <typename LineStringIterator>
double rangedLength(LineStringIterator start, LineStringIterator end) {
  double l = 0.;
  helper::forEachPair(start, end, [&l](const auto& seg1, const auto& seg2) { l += distance(seg1, seg2); });
  return l;
}

template <typename LineStringT>
std::vector<double> lengthRatios(const LineStringT& lineString) {
  std::vector<double> lengths;
  if (lineString.size() <= 1) {
    return lengths;
  }
  if (lineString.size() == 2) {
    return {1.};
  }
  const auto totalLength = length(lineString);
  lengths.reserve(lineString.size() - 1);
  helper::forEachPair(lineString.begin(), lineString.end(), [&lengths, &totalLength](const auto& p1, const auto& p2) {
    lengths.push_back(distance(p1, p2) / totalLength);
  });
  return lengths;
}

template <typename LineStringT>
std::vector<double> accumulatedLengthRatios(const LineStringT& lineString) {
  auto lengths = lengthRatios(lineString);
  helper::forEachPair(lengths.begin(), lengths.end(), [](const auto& l1, auto& l2) { l2 += l1; });
  return lengths;
}

template <typename LineStringT>
traits::BasicPointT<traits::PointType<LineStringT>> interpolatedPointAtDistance(LineStringT lineString, double dist) {
  assert(!lineString.empty());
  if (dist < 0) {
    lineString = internal::invert(lineString);
    dist = -dist;
  }

  double currentCumulativeLength = 0.0;
  for (auto first = lineString.begin(), second = std::next(lineString.begin()); second != lineString.end();
       ++first, ++second) {
    const auto p1 = traits::toBasicPoint(*first);
    const auto p2 = traits::toBasicPoint(*second);
    double currentLength = distance(p1, p2);
    currentCumulativeLength += currentLength;
    if (currentCumulativeLength >= dist) {
      double remainingDistance = dist - (currentCumulativeLength - currentLength);
      if (remainingDistance < 1.e-8) {
        return p1;
      }
      return p1 + remainingDistance / currentLength * (p2 - p1);
    }
  }
  return traits::toBasicPoint(lineString.back());
}

template <typename LineStringT>
traits::PointType<LineStringT> nearestPointAtDistance(LineStringT lineString, double dist) {
  using traits::toBasicPoint;
  assert(!lineString.empty());
  if (dist < 0) {
    lineString = internal::invert(lineString);
    dist = -dist;
  }
  double currentCumulativeLength = 0.0;
  for (auto first = lineString.begin(), second = std::next(lineString.begin()); second != lineString.end();
       ++first, ++second) {
    const auto& p1 = *first;
    const auto& p2 = *second;
    double currentLength = distance(p1, p2);
    currentCumulativeLength += currentLength;
    if (currentCumulativeLength >= dist) {
      double remainingDistance = dist - (currentCumulativeLength - currentLength);
      if (remainingDistance > currentLength / 2.0) {
        return p2;
      }
      return p1;
    }
  }
  return lineString.back();
}

template <typename LineString3dT>
double signedDistance(const LineString3dT& lineString, const BasicPoint3d& p) {
  static_assert(traits::is3D<LineString3dT>(), "Please call this function with a 3D type!");
  return internal::signedDistanceImpl(lineString, p).first;
}

template <typename LineString2dT>
double signedDistance(const LineString2dT& lineString, const BasicPoint2d& p) {
  static_assert(traits::is2D<LineString2dT>(), "Please call this function with a 2D type!");
  return internal::signedDistanceImpl(lineString, p).first;
}

template <typename LineString2dT>
ArcCoordinates toArcCoordinates(const LineString2dT& lineString, const BasicPoint2d& point) {
  auto res = internal::signedDistanceImpl(lineString, point);
  auto dist = res.first;
  const auto& projectedPoint = res.second;
  // find first point in segment in linestring
  double length = 0.;
  auto accumulateLength = [&length, &point = projectedPoint.result->segmentPoint1](const auto& first,
                                                                                   const auto& second) {
    if (boost::geometry::equals(first, point)) {
      return true;
    }
    length += distance(first, second);
    return false;
  };
  helper::forEachPairUntil(lineString.begin(), lineString.end(), accumulateLength);
  length += distance(projectedPoint.result->segmentPoint1, projectedPoint.result->projectedPoint);
  return {length, dist};
}

template <typename LineString3dT>
IfLS<LineString3dT, BoundingBox3d> boundingBox3d(const LineString3dT& lineString) {
  static_assert(traits::is3D<LineString3dT>(), "Please call this function with a 3D type!");
  BoundingBox3d bb;
  for (const auto& p : lineString) {
    bb.extend(traits::toBasicPoint(p));
  }
  return bb;
}

template <typename LineString2dT>
IfLS<LineString2dT, BoundingBox2d> boundingBox2d(const LineString2dT& lineString) {
  BoundingBox2d bb;
  for (const auto& p : traits::to2D(lineString)) {
    bb.extend(traits::toBasicPoint(p));
  }
  return bb;
}

template <typename LineString3dT, typename>
BasicPoint3d project(const LineString3dT& lineString, const BasicPoint3d& pointToProject) {
  static_assert(traits::is3D<LineString3dT>(), "Please call this function with a 3D type!");
  helper::ProjectedPoint<BasicPoint3d> projectedPoint;
  distance(lineString, pointToProject, projectedPoint);
  return projectedPoint.result->projectedPoint;
}

template <typename LineString2dT, typename>
BasicPoint2d project(const LineString2dT& lineString, const BasicPoint2d& pointToProject) {
  static_assert(traits::is2D<LineString2dT>(), "Please call this function with a 2D type!");
  helper::ProjectedPoint<BasicPoint2d> projectedPoint;
  distance(lineString, pointToProject, projectedPoint);
  return projectedPoint.result->projectedPoint;
}

template <typename LineString3dT>
IfLS<LineString3dT, bool> intersects3d(const LineString3dT& linestring, const LineString3dT& otherLinestring,
                                       double heightTolerance) {
  auto ls2d(traits::toHybrid(traits::to2D(linestring)));
  auto ols2d(traits::toHybrid(traits::to2D(otherLinestring)));
  BasicPoints2d intersections;
  boost::geometry::intersection(ls2d, ols2d, intersections);
  auto distanceSmallerTolerance = [heightTolerance, &linestring, &otherLinestring](const auto& point) {
    auto pProj1 = project(linestring, BasicPoint3d(point.x(), point.y(), 0));
    auto pProj2 = project(otherLinestring, BasicPoint3d(point.x(), point.y(), 0));
    return distance(pProj1, pProj2) < heightTolerance;
  };
  return std::any_of(intersections.begin(), intersections.end(), distanceSmallerTolerance);
}

template <typename LineString3dT>
IfLS<LineString3dT, std::pair<BasicPoint3d, BasicPoint3d>> projectedPoint3d(const LineString3dT& l1,
                                                                            const LineString3dT& l2) {
  return internal::projectedPoint3d(traits::toHybrid(l1), traits::toHybrid(l2));
}

template <typename LineStringT>
IfLS<LineStringT, double> distance2d(const LineStringT& l1, const LineStringT& l2) {
  return distance(traits::toHybrid(traits::to2D(l1)), traits::toHybrid(traits::to2D(l2)));
}

template <typename LineStringT>
IfLS<LineStringT, double> distance2d(const LineStringT& l1, const BasicPoint2d& p) {
  return distance(traits::to2D(traits::toConst(l1)), p);
}

template <typename LineString3dT>
IfLS<LineString3dT, double> distance3d(const LineString3dT& l1, const LineString3dT& l2) {
  auto projPoint = internal::projectedPoint3d(traits::toHybrid(l1), traits::toHybrid(l2));
  return (projPoint.first - projPoint.second).norm();
}

template <typename LineString1T, typename LineString2T>
std::pair<LineString1T, LineString2T> align(LineString1T left, LineString2T right) {
  using traits::toBasicPoint;
  // degenerated case
  if ((left.size() <= 1 && right.size() <= 1) || right.empty() || left.empty()) {
    return {left, right};
  }
  auto getMiddlePoint = [](auto& ls) {
    return ls.size() > 2 ? toBasicPoint(ls[ls.size() / 2]) : (toBasicPoint(ls.front()) + toBasicPoint(ls.back())) * 0.5;
  };
  //! @todo this sadly is a bit heuristical...
  bool rightOfLeft = signedDistance(left, getMiddlePoint(right)) < 0;
  if (!rightOfLeft && left.size() > 1) {
    left = left.invert();
  }

  bool leftOfRight = signedDistance(right, getMiddlePoint(left)) > 0;
  if (!leftOfRight && right.size() > 1) {
    right = right.invert();
  }
  return {left, right};
}

template <typename LineString2dT>
BasicPoint2d fromArcCoordinates(const LineString2dT& lineString, const ArcCoordinates& arcCoords) {
  auto hLineString = utils::toHybrid(lineString);
  auto ratios = accumulatedLengthRatios(lineString);
  const auto llength = length(lineString);
  int startIdx = -1, endIdx = -1;
  for (int i = 0; i < static_cast<int>(ratios.size()); ++i) {
    if (ratios.at(static_cast<size_t>(i)) * llength > arcCoords.length) {
      startIdx = i;
      endIdx = i + 1;
      break;
    }
  }
  if (endIdx == -1) {
    endIdx = lineString.size() - 1;
    startIdx = endIdx - 1;
  }
  return internal::fromArcCoords(hLineString, interpolatedPointAtDistance(utils::to2D(lineString), arcCoords.length),
                                 startIdx, endIdx, arcCoords.distance);
}

// negative means from back
template <typename LineString2dT>
BasicPoint2d fromArcCoordinatesAtPoint(const LineString2dT& lineString, const int idx, const double distance) {
  if (std::abs(idx) > lineString.size() - 1) throw InvalidInputError("Index out of bounds");
  int startIdx = (idx >= 0) ? std::max(0, idx - 1) : std::max(0, static_cast<int>(lineString.size()) + idx - 1);
  int endIdx = (idx >= 0)
                   ? std::min(idx + 1, static_cast<int>(lineString.size()) - 1)
                   : std::min(static_cast<int>(lineString.size()) + idx + 1, static_cast<int>(lineString.size()) - 1);
  int pIdx = (idx >= 0) ? idx : static_cast<int>(lineString.size()) + idx;
  auto hLineString = utils::toHybrid(lineString);
  return internal::fromArcCoords(hLineString, hLineString[pIdx], startIdx, endIdx, distance);
}

// template <typename LineString2dT>
// SelfIntersections2d getSelfIntersections(const LineString2dT& lineString) {
//  if (lineString.size() <= 3) return SelfIntersections2d();
//  SelfIntersections2d ret;
//  auto lsFromSegment = [](const typename LineString2dT::ConstSegmentType& seg) -> BasicLineString2d {
//    return BasicLineString2d({utils::toBasicPoint(seg.first), utils::toBasicPoint(seg.second)});
//  };
//  auto intersectSegments = [lsFromSegment](
//                               const typename LineString2dT::ConstSegmentType& lhs,
//                               const typename LineString2dT::ConstSegmentType& rhs) -> boost::optional<BasicPoint2d> {
//    auto ls1 = lsFromSegment(lhs);
//    auto ls2 = lsFromSegment(rhs);
//    boost::optional<BasicPoint2d> res;
//    BasicPoints2d intersection;
//    boost::geometry::intersection(ls1, ls2, intersection);
//    assert(intersection.size() < 2);
//    if (intersection.empty()) return res;
//    res = intersection.front();
//    return res;
//  };
//  for (int i = 0; i < lineString.size() - 1; ++i) {
//    const auto& segment = lineString.segment(i);
//    for (int j = i + 2; j < lineString.size() - 1; ++j) {
//      const auto& otherSegment = lineString.segment(j);
//      auto intersection = intersectSegments(segment, otherSegment);
//      if (intersection) {
//        ret.emplace_back(SelfIntersection2d{i, j, *intersection});
//      }
//    }
//  }
//  std::sort(
//      ret.begin(), ret.end(), [&lineString](const SelfIntersection2d& lhs, const SelfIntersection2d& rhs) -> bool {
//        if (lhs.firstSegmentIdx == rhs.firstSegmentIdx) {
//          BasicPoints2d pointsToCompare({lhs.intersectionPoint, rhs.intersectionPoint});
//          if (internal::sortAlongSIdxs(BasicLineString2d({utils::toBasicPoint(lineString[lhs.firstSegmentIdx]),
//                                                          utils::toBasicPoint(lineString[lhs.firstSegmentIdx + 1])}),
//                                       pointsToCompare)
//                  .front() == 0)
//            return true;
//          return false;
//        }
//        return (lhs.firstSegmentIdx < rhs.firstSegmentIdx);
//      });
//  return ret;
//}

template <typename LineString2dT>
SelfIntersections2d getSelfIntersections(const LineString2dT& lineString) {
  if (lineString.size() <= 3) return SelfIntersections2d();
  SelfIntersections2d ret;
  using Segment = std::pair<size_t, size_t>;
  auto intersectSegments = [&lineString](const Segment& lhs, const Segment& rhs) -> boost::optional<BasicPoint2d> {
    auto ls1 = BasicLineString2d({lineString.at(lhs.first), lineString.at(lhs.second)});
    auto ls2 = BasicLineString2d({lineString.at(rhs.first), lineString.at(rhs.second)});
    boost::optional<BasicPoint2d> res;
    BasicPoints2d intersection;
    boost::geometry::intersection(ls1, ls2, intersection);
    assert(intersection.size() < 2);
    if (intersection.empty()) return res;
    res = intersection.front();
    return res;
  };
  for (size_t i = 0; i < lineString.size() - 1; ++i) {
    const auto& segment = std::make_pair(i, i + 1);
    for (size_t j = i + 2; j < lineString.size() - 1; ++j) {
      const auto& otherSegment = std::make_pair(j, j + 1);
      auto intersection = intersectSegments(segment, otherSegment);
      if (intersection) {
        ret.emplace_back(SelfIntersection2d{static_cast<int>(i), static_cast<int>(j), *intersection});
      }
    }
  }
  std::sort(
      ret.begin(), ret.end(), [&lineString](const SelfIntersection2d& lhs, const SelfIntersection2d& rhs) -> bool {
        if (lhs.firstSegmentIdx == rhs.firstSegmentIdx) {
          BasicPoints2d pointsToCompare({lhs.intersectionPoint, rhs.intersectionPoint});
          if (internal::sortAlongSIdxs(BasicLineString2d({utils::toBasicPoint(lineString[lhs.firstSegmentIdx]),
                                                          utils::toBasicPoint(lineString[lhs.firstSegmentIdx + 1])}),
                                       pointsToCompare)
                  .front() == 0)
            return true;
          return false;
        }
        return (lhs.firstSegmentIdx < rhs.firstSegmentIdx);
      });
  return ret;
}

template <typename LineString2dT>
BasicLineString2d offset(const LineString2dT& lineString, const double distance) {
  if (lineString.size() < 1) throw std::runtime_error("Not shifting linestring with less than two points");
  BasicLineString2d ret(lineString.size());
  for (int i = 0; i < lineString.size(); ++i) {
    ret[i] = internal::shiftLateral(lineString, i, distance);
  }

  auto earliestIntersectionBehind = [&lineString](const SelfIntersections2d& sInt, const int intIdx) {
    using SPoint = std::tuple<int, double>;
    using SList = std::vector<SPoint>;
    const auto segToCheckIdx = sInt.at(intIdx).lastSegmentIdx;
    BasicLineString2d segLs(
        {utils::toBasicPoint(lineString[segToCheckIdx]), utils::toBasicPoint(lineString[segToCheckIdx + 1])});
    SList startingIntersections;
    for (int i = 0; i < sInt.size(); ++i) {
      if (sInt.at(i).firstSegmentIdx == segToCheckIdx)
        startingIntersections.emplace_back(
            std::make_tuple(i, toArcCoordinates(segLs, sInt.at(i).intersectionPoint).length));
    }
    std::sort(startingIntersections.begin(), startingIntersections.end(),
              [](const SPoint& lhs, const SPoint& rhs) { return std::get<1>(lhs) < std::get<1>(rhs); });
    auto sToBeat = toArcCoordinates(segLs, sInt.at(intIdx).intersectionPoint).length;
    for (const auto& si : startingIntersections) {
      if (std::get<1>(si) > sToBeat) return std::get<0>(si);
    }
    return -1;
  };
  auto selfIntersections = getSelfIntersections(ret);
  if (selfIntersections.empty()) return ret;
  int curSegment = 0;
  int nextIntersection = 0;
  BasicLineString2d res({ret[0]});
  while (nextIntersection < selfIntersections.size()) {
    res.insert(res.end(), ret.begin() + curSegment,
               ret.begin() + selfIntersections.at(nextIntersection).firstSegmentIdx);
    res.push_back(selfIntersections.at(nextIntersection).intersectionPoint);
    curSegment = selfIntersections.at(nextIntersection).lastSegmentIdx;
    auto earliestOnSameSegment = earliestIntersectionBehind(selfIntersections, nextIntersection);
    if (earliestOnSameSegment != -1) {
      nextIntersection = earliestOnSameSegment;
    }
    while (nextIntersection < selfIntersections.size() &&
           selfIntersections.at(nextIntersection).firstSegmentIdx < curSegment)
      ++nextIntersection;
  }
  if (curSegment < ret.size() - 1) res.insert(res.end(), ret.begin() + curSegment + 1, ret.end());
  return res;
}
}  // namespace geometry

}  // namespace lanelet
