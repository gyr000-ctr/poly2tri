/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <list>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "../common/shapes.h"

namespace p2t {

// Inital triangle factor, seed triangle will extend 30% of
// PointSet width to both left and right.
const double kAlpha = 0.3;

struct Point;
class Triangle;
struct Node;
struct Edge;
class AdvancingFront;

/// Classification of an edge in the produced triangulation.
enum EdgeType {
  /// Edge of the outer polyline (the original outline)
  EDGE_OUTLINE = 0,
  /// Edge of a hole polyline
  EDGE_HOLE = 1,
  /// Any other edge (created by the triangulation, including edges
  /// incident to interior Steiner points)
  EDGE_INTERIOR = 2
};

/// A unique edge of the triangulation together with the point indices of
/// its two endpoints. Indices refer to the internal point list, i.e. the
/// order returned by SweepContext::GetPoint() after Triangulate() (sorted
/// by ascending y, then ascending x). The same pair of points always maps
/// to the same two indices.
struct EdgeInfo {
  /// Index of the first endpoint (lower y, then lower x)
  std::size_t p;
  /// Index of the second endpoint
  std::size_t q;
  /// Outline / hole / interior classification
  EdgeType type;
};

namespace internal {

// Canonical, order-independent key of an undirected point pair.
typedef std::pair<const Point*, const Point*> PointPairKey;

struct PointPairHash {
  std::size_t operator()(const PointPairKey& key) const
  {
    // Order the two pointers so (a, b) and (b, a) hash identically.
    const Point* lo = key.first;
    const Point* hi = key.second;
    if (lo > hi) {
      std::swap(lo, hi);
    }
    const std::size_t h1 = std::hash<const Point*>()(lo);
    const std::size_t h2 = std::hash<const Point*>()(hi);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
  }
};

struct PointPairEqual {
  bool operator()(const PointPairKey& a, const PointPairKey& b) const
  {
    const Point* a_lo = a.first;
    const Point* a_hi = a.second;
    if (a_lo > a_hi) {
      std::swap(a_lo, a_hi);
    }
    const Point* b_lo = b.first;
    const Point* b_hi = b.second;
    if (b_lo > b_hi) {
      std::swap(b_lo, b_hi);
    }
    return a_lo == b_lo && a_hi == b_hi;
  }
};

} // namespace internal

class SweepContext {
public:

/// Constructor
explicit SweepContext(std::vector<Point*> polyline);
/// Destructor
~SweepContext();

void set_head(Point* p1);

Point* head() const;

void set_tail(Point* p1);

Point* tail() const;

size_t point_count() const;

Node* LocateNode(const Point& point);

void RemoveNode(Node* node);

void CreateAdvancingFront();

/// Try to map a node to all sides of this triangle that don't have a neighbor
void MapTriangleToNodes(Triangle& t);

void AddToMap(Triangle* triangle);

Point* GetPoint(size_t index);

Point* GetPoints();

void RemoveFromMap(Triangle* triangle);

void AddHole(const std::vector<Point*>& polyline);

void AddPoint(Point* point);

AdvancingFront* front() const;

void MeshClean(Triangle& triangle);

std::vector<Triangle*> &GetTriangles();
std::list<Triangle*> &GetMap();

/// Number of points used by the triangulation (polyline points, hole
/// points and Steiner points). Valid after Triangulate().
std::size_t GetPointCount() const;

/// All unique edges of the triangulation, each classified as an outline
/// edge (outer polyline), a hole edge or an interior edge, with the point
/// indices of its endpoints. Sorted by ascending (p, q) index pairs.
/// Valid after Triangulate().
const std::vector<EdgeInfo>& GetEdges();

/// Only the outline edges, i.e. the edges of the original outer polyline,
/// with the point indices of their endpoints. Sorted by ascending
/// (p, q) index pairs. Valid after Triangulate().
std::vector<EdgeInfo> GetOutlineEdges();

std::vector<Edge*> edge_list;

struct Basin {
  Node* left_node;
  Node* bottom_node;
  Node* right_node;
  double width;
  bool left_highest;

  Basin()
  : left_node(nullptr), bottom_node(nullptr), right_node(nullptr), width(0.0), left_highest(false)
  {
  }

  void Clear()
  {
    left_node = nullptr;
    bottom_node = nullptr;
    right_node = nullptr;
    width = 0.0;
    left_highest = false;
  }
};

struct EdgeEvent {
  Edge* constrained_edge;
  bool right;

  EdgeEvent() : constrained_edge(NULL), right(false)
  {
  }
};

Basin basin;
EdgeEvent edge_event;

private:

friend class Sweep;

std::vector<Triangle*> triangles_;
std::list<Triangle*> map_;
std::vector<Point*> points_;
std::unordered_set<Point*> point_set_;

// Advancing front
AdvancingFront* front_;
// head point used with advancing front
Point* head_;
// tail point used with advancing front
Point* tail_;

Node *af_head_, *af_middle_, *af_tail_;

// Constraint edges belonging to hole polylines, used to distinguish hole
// edges from outer polyline edges while classifying the triangulation.
std::unordered_set<internal::PointPairKey, internal::PointPairHash,
  internal::PointPairEqual> hole_edges_;

// Lazily built edge classification, rebuilt on every Triangulate().
std::vector<EdgeInfo> edge_infos_;
bool edge_infos_dirty_;

void InitTriangulation();
void InitEdges(const std::vector<Point*>& polyline);

};

inline AdvancingFront* SweepContext::front() const
{
  return front_;
}

inline size_t SweepContext::point_count() const
{
  return points_.size();
}

inline void SweepContext::set_head(Point* p1)
{
  head_ = p1;
}

inline Point* SweepContext::head() const
{
  return head_;
}

inline void SweepContext::set_tail(Point* p1)
{
  tail_ = p1;
}

inline Point* SweepContext::tail() const
{
  return tail_;
}

inline std::size_t SweepContext::GetPointCount() const
{
  return points_.size();
}

inline const std::vector<EdgeInfo>& SweepContext::GetEdges()
{
  if (edge_infos_dirty_) {
    // Map every point pointer to its stable index in the internal point
    // list (ascending y, then ascending x after the sweep's sort).
    std::unordered_map<const Point*, std::size_t> point_index;
    point_index.reserve(points_.size() * 2);
    for (std::size_t i = 0; i < points_.size(); ++i) {
      point_index[points_[i]] = i;
    }

    // Every undirected edge is visited once per incident triangle; the
    // set collapses interior edges (two triangles) to a single entry.
    std::unordered_set<internal::PointPairKey, internal::PointPairHash,
      internal::PointPairEqual> seen;
    seen.reserve(triangles_.size() * 3);

    edge_infos_.clear();
    for (Triangle* t : triangles_) {
      for (int i = 0; i < 3; ++i) {
        const Point* a = t->GetPoint(i);
        const Point* b = t->GetPoint((i + 1) % 3);
        internal::PointPairKey key(a, b);
        if (!seen.insert(key).second) {
          continue;
        }
        // constrained_edge[k] belongs to the edge opposite points_[k],
        // i.e. the edge (points_[(k+1)%3], points_[(k+2)%3]); for the
        // edge (points_[i], points_[(i+1)%3]) that index is (i+2)%3.
        const int opposite = (i + 2) % 3;
        EdgeType type;
        if (t->constrained_edge[opposite]) {
          // A constrained edge is either part of a hole polyline or of
          // the outer polyline (the original outline).
          type = hole_edges_.count(key) ? EDGE_HOLE : EDGE_OUTLINE;
        } else {
          type = EDGE_INTERIOR;
        }
        std::size_t ia = point_index[a];
        std::size_t ib = point_index[b];
        EdgeInfo info;
        info.p = ia < ib ? ia : ib;
        info.q = ia < ib ? ib : ia;
        info.type = type;
        edge_infos_.push_back(info);
      }
    }

    std::sort(edge_infos_.begin(), edge_infos_.end(),
      [](const EdgeInfo& a, const EdgeInfo& b) {
        if (a.p != b.p) {
          return a.p < b.p;
        }
        return a.q < b.q;
      });
    edge_infos_dirty_ = false;
  }
  return edge_infos_;
}

inline std::vector<EdgeInfo> SweepContext::GetOutlineEdges()
{
  std::vector<EdgeInfo> outline;
  for (const EdgeInfo& edge : GetEdges()) {
    if (edge.type == EDGE_OUTLINE) {
      outline.push_back(edge);
    }
  }
  return outline;
}

}
