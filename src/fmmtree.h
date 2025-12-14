#pragma once

#include "point.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

enum class NodeType : uint8_t { Internal, Leaf };

template <class Kernel> struct FmmNode2 {
  using Multipole = typename Kernel::Multipole;
  using Local = typename Kernel::Local;

  FmmNode2(int p, Box2 box, NodeType type)
      : box(box), type(type), multipole(p, box.center), local(p, box.center) {}

  Box2 box;
  NodeType type;
  std::vector<Point> points; // for leaf nodes

  FmmNode2 *parent = nullptr;
  FmmNode2 *children[4] = {nullptr, nullptr, nullptr, nullptr};

  Multipole multipole;
  Local local;

  std::vector<FmmNode2 *> near_neighbors;
  std::vector<FmmNode2 *> interaction_list;

  bool is_leaf() const { return type == NodeType::Leaf; }
};

// Balanced FMM tree implementation
template <class Kernel> struct NaiveFmmTree {
  using Multipole = typename Kernel::Multipole;
  using Local = typename Kernel::Local;

  using Node = FmmNode2<Kernel>;

  NaiveFmmTree(int p, const std::vector<Point> &sources, int height);
  ~NaiveFmmTree();

  Node *root() const { return root_; }
  int height() const { return height_; }
  std::size_t num_leaves() const { return std::pow(4, height_); }

  double evaluate(Vector2 point) const;
  std::vector<double> evaluateSources() const;

protected:
  int p_;
  Node *root_;
  int height_;
  std::vector<std::vector<Node *>> levels_;
  const std::vector<Point> &sources_;

private:
  std::size_t getLeafIndex(Vector2 position) const;
  void buildChildNodes(Node *node, int level);
  void computeNodeLists(Node *node);
};

template <class Kernel>
std::size_t NaiveFmmTree<Kernel>::getLeafIndex(Vector2 position) const {
  std::size_t index = 0;
  Box2 bbox = root_->box;
  for (int level = 0; level < height_; level++) {
    uint32_t quadrant = getQuadrant(bbox, position);
    index = (index << 2) | quadrant;
    bbox = getChildBox(bbox, quadrant);
  }
  return index;
}

template <class Kernel>
NaiveFmmTree<Kernel>::NaiveFmmTree(int p, const std::vector<Point> &sources,
                                   int height)
    : p_(p), height_(height), sources_(sources) {
  if (height_ == 0) {
    return;
  }
  levels_.resize(height_ + 1);

  Box2 root_box = computeBoundingBox(sources);
  root_ = new Node(p_, root_box, NodeType::Leaf);
  levels_[0].push_back(root_);

  for (int level = 0; level <= height_; level++) {
    for (Node *node : levels_[level]) {
      buildChildNodes(node, level);
      computeNodeLists(node);
    }
  }

  for (const Point &p : sources) {
    std::size_t leafIndex = getLeafIndex(p.position);
    levels_[height_][leafIndex]->points.push_back(p);
  }

  // step 1: form multipole expansions at each leaf node
  for (std::size_t i = 0; i < num_leaves(); i++) {
    Node *leaf = levels_[height_][i];
    leaf->multipole.buildExpansion(leaf->points);
  }

  // step 2: form multipole expansions up the tree by combining child multipole
  // expansions
  for (int level = height_ - 1; level >= 2; level--) {
    for (Node *node : levels_[level]) {
      Multipole me(p_, node->box.center);

      for (Node *child : node->children) {
        Vector2 shift = child->box.center - node->box.center;
        me += child->multipole.M2M(Complex(shift.x, shift.y));
      }
      node->multipole = me;
    }
  }

  // step 3: form local expansions down the tree by combining multipole
  // expansions from interaction list
  for (int level = 2; level <= height_; level++) {
    for (Node *node : levels_[level]) {
      Vector2 shift = node->parent->box.center - node->box.center;
      node->local += node->parent->local.L2L(Complex(shift.x, shift.y));

      for (Node *interaction : node->interaction_list) {
        Local le(p_, node->box.center);
        le.M2L(interaction->multipole);
        node->local += le;
      }
    }
  }
}

template <class Kernel> NaiveFmmTree<Kernel>::~NaiveFmmTree() {
  for (int level = 0; level <= height_; level++) {
    for (Node *node : levels_[level]) {
      delete node;
    }
  }
}

template <class Kernel>
void NaiveFmmTree<Kernel>::buildChildNodes(Node *node, int level) {
  if (level == height_) {
    node->type = NodeType::Leaf;
    return;
  }

  node->type = NodeType::Internal;
  for (int q = 0; q < 4; q++) {
    Box2 child_box = getChildBox(node->box, q);
    Node *child = new Node(p_, child_box, NodeType::Leaf);
    child->parent = node;
    node->children[q] = child;
    levels_[level + 1].push_back(child);
  }
}

// compute near neighbors and interaction list for node, assuming parent node
// has already been computed
template <class Kernel>
void NaiveFmmTree<Kernel>::computeNodeLists(Node *node) {
  Node *parent = node->parent;
  if (parent == nullptr) {
    node->near_neighbors.push_back(node);
    return;
  }

  for (Node *parent_neighbor : parent->near_neighbors) {
    for (Node *child : parent_neighbor->children) {
      if (adjacent(node->box, child->box)) {
        node->near_neighbors.push_back(child);
      } else {
        node->interaction_list.push_back(child);
      }
    }
  }
}

template <class Kernel>
double NaiveFmmTree<Kernel>::evaluate(Vector2 point) const {
  std::size_t leaf_index = getLeafIndex(point);
  Node *leaf = levels_[height_][leaf_index];

  double result = leaf->local.evaluate(point);

  for (Node *near_neighbor : leaf->near_neighbors) {
    for (const Point &p : near_neighbor->points) {
      if (p.position == point) {
        continue;
      }
      result += Kernel::potential(p, point);
    }
  }

  return result;
}

template <class Kernel>
std::vector<double> NaiveFmmTree<Kernel>::evaluateSources() const {
  int num_sources = sources_.size();
  std::vector<double> potentials(num_sources);

  for (int i = 0; i < num_sources; i++) {
    potentials[i] = evaluate(sources_[i].position);
  }

  return potentials;
}
