// Copysecond

#ifndef SOURCE_RTREE_HPP_
#define SOURCE_RTREE_HPP_

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <queue>
#include "Rectangle.hpp"

template <size_t N, typename ElemType, size_t M, size_t m = M / 2>
class RTree {
 public:
  struct Node;

  struct SpatialObject {
    Rectangle<N> box;
    ElemType identifier;
    std::shared_ptr<Node> child_pointer;
  };

  struct Node {
    typedef SpatialObject *iterator;
    typedef const SpatialObject *const_iterator;

    iterator begin();
    iterator end();

    const_iterator begin() const;
    const_iterator end() const;

    SpatialObject &operator[](size_t index);
    SpatialObject operator[](size_t index) const;

    bool is_leaf();

    std::shared_ptr<Node> insert(const SpatialObject &new_entry);
    void pickSeeds(std::vector<SpatialObject> &list, int &a, int &b);
    int pickNext(std::vector<SpatialObject> &input_1,
                 Rectangle<N> &bb_first,
                 Rectangle<N> &bb_second);
    
    SpatialObject entry[M];
    size_t size = 0;
  };

  RTree();
  virtual ~RTree();
  size_t dimension() const;
  size_t size() const;
  bool empty() const;
  void print();

  void insert(const Rectangle<N> &box, const ElemType &value);
  std::shared_ptr<Node> choose_leaf(const std::shared_ptr<Node> &current_node,
                                    const Rectangle<N> &box,
                                    const ElemType &value);

  std::shared_ptr<Node> choose_node(const std::shared_ptr<Node> &current_node,
                                    const Rectangle<N> &box,
                                    SpatialObject *&entry);

  std::shared_ptr<Node> adjust_tree(const std::shared_ptr<Node> &parent,
                                    const std::shared_ptr<Node> &first,
                                    const std::shared_ptr<Node> &second,
                                    SpatialObject *entry);

  std::vector<ElemType> operator[](const Rectangle<N> &box);
  std::vector<ElemType> at(const Rectangle<N> &box);
  const std::vector<ElemType> at(const Rectangle<N> &box) const;
  // std::vector<ElemType> kNNValue(const Rectangle<N> &box, size_t k) const;

  void search(const std::shared_ptr<Node> &cur_node,
              const Rectangle<N> &box,
              std::vector<ElemType> &value);

  //private:
  std::shared_ptr<Node> root_pointer_;
  size_t entries;
};

/** Node R-tree struct implementation details*/

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::iterator
RTree<N, ElemType, M, m>::Node::begin() {
  return entry;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::iterator
RTree<N, ElemType, M, m>::Node::end() {
  return entry + size;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::const_iterator
RTree<N, ElemType, M, m>::Node::begin() const {
  return entry;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::Node::const_iterator
RTree<N, ElemType, M, m>::Node::end() const {
  return entry + size;
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::SpatialObject
    &RTree<N, ElemType, M, m>::Node::operator[](size_t index) {
  return entry[index];
}

template <size_t N, typename ElemType, size_t M, size_t m>
typename RTree<N, ElemType, M, m>::SpatialObject
RTree<N, ElemType, M, m>::Node::operator[](size_t index) const {
  return entry[index];
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType> RTree<N, ElemType, M, m>::operator[](const Rectangle<N> &box) {
  std::vector<ElemType> value;
  search(root_pointer_, box, value);
  return value;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::vector<ElemType> RTree<N, ElemType, M, m>::at(const Rectangle<N> &box){
  std::vector<ElemType> value;
  search(root_pointer_, box, value);
  return value;
}

template <size_t N, typename ElemType, size_t M, size_t m>
const std::vector<ElemType> RTree<N, ElemType, M, m>::at(const Rectangle<N> &box) const {
  std::vector<ElemType> value;
  search(root_pointer_, box, value);
  return value;
}
template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::Node::is_leaf() {
  if (size && entry[0].child_pointer) {
    return false;
  }
  return true;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::Node::insert(const SpatialObject &new_entry) {
  if (size < M) {
    entry[size++] = new_entry;
    return nullptr;
  }
  
 
  // TODO(ADE): Split the entries and return a pointer to new node
  // caused due to split.
  std::shared_ptr<Node> new_node = std::make_shared<Node>();
  std::vector<SpatialObject> split_list(M + 1);
  for (size_t i = 0; i < M; ++i) {
    split_list[i] = entry[i];
  }
  split_list[M] = new_entry;
  this->size = 0;
  
  int first, second;
  pickSeeds(split_list, first, second);
  this->insert(split_list[first]);
  new_node->insert(split_list[second]);
  if(first > second) std::swap(first, second);
  split_list.erase(split_list.begin() + first);
  split_list.erase(split_list.begin() + second - 1);
  
  Rectangle<N> bb_first = (*this)[0].box;
  Rectangle<N> bb_second = (*new_node)[0].box;
  while (!split_list.empty()) {
    if (this->size + split_list.size() == m) {
      for (SpatialObject &obj : split_list) {
        this->insert(obj);
      }
      split_list.clear();
      continue;
    }
    if (new_node->size + split_list.size() == m) {
      for (SpatialObject &obj : split_list) {
        new_node->insert(obj);
      }
      split_list.clear();
      continue;
    }
    int next = pickNext(split_list, bb_first, bb_second);
    Rectangle<N> enlargement_first = split_list[next].box;
    Rectangle<N> enlargement_second = split_list[next].box;
    enlargement_first.adjust(bb_first);
    enlargement_second.adjust(bb_second);
    float area_first = enlargement_first.get_area() - split_list[next].box.get_area();
    float area_second = enlargement_second.get_area() - split_list[next].box.get_area();
    if (area_first < area_second ||
        (area_first == area_second && bb_first.get_area() < bb_second.get_area()) ||
        (area_first == area_second && bb_first.get_area() == bb_second.get_area() &&
         this->size < new_node->size)) {
      bb_first.adjust(split_list[next].box);
      this->insert(split_list[next]);
    }
    else {
      bb_second.adjust(split_list[next].box);
      new_node->insert(split_list[next]);
    }
    split_list.erase(split_list.begin() + next);
  }
  return new_node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::Node::pickSeeds(std::vector<SpatialObject> &input_1,int &seed_1, int &seed_2) {
  seed_1 = 0, seed_2 = 1;
  float area;
  Rectangle<N> temp_box = input_1[seed_1].box;
  float limit_area = temp_box.get_area();
  temp_box.adjust(input_1[seed_2].box);
  limit_area -= (input_1[seed_1].box.get_area() + input_1[seed_2].box.get_area());
  for (int i = 0; i < M; ++i) {
    for (int j = i+1; j <= M; ++j) {
      temp_box = input_1[i].box;
      temp_box.adjust(input_1[j].box);
      area = temp_box.get_area() - input_1[i].box.get_area() - input_1[j].box.get_area();
      if (area > limit_area) {
        limit_area = area;
        seed_1 = i;
        seed_2 = j;
      }
    }
  }
}

template <size_t N, typename ElemType, size_t M, size_t m>
int RTree<N, ElemType, M, m>::Node::pickNext(std::vector<SpatialObject> &input_1,
                                             Rectangle<N> &bb_first,
                                             Rectangle<N> &bb_second) {
  int ind = 0;
  Rectangle<N> enlargement_first = input_1[0].box;
  Rectangle<N> enlargement_second = input_1[0].box;
  enlargement_first.adjust(bb_first);
  enlargement_second.adjust(bb_second);
  float area_first = enlargement_first.get_area() - input_1[0].box.get_area();
  float area_second = enlargement_second.get_area() - input_1[0].box.get_area();
  float max_d = abs(area_first - area_second);
  float area_d;
  for (int i = 1; i < input_1.size(); ++i) {
    enlargement_first = input_1[i].box;
    enlargement_second = input_1[i].box;
    enlargement_first.adjust(bb_first);
    enlargement_second.adjust(bb_second);
    area_first = enlargement_first.get_area() - input_1[i].box.get_area();
    area_second = enlargement_second.get_area() - input_1[i].box.get_area();
    area_d = abs(area_first - area_second);
    if (area_d > max_d) {
      max_d = area_d;
      ind = i;
    }
  }
  return ind;
}

/** R-Tree class implementation details */

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::RTree() : root_pointer_(std::make_shared<Node>()), entries(0) {}

template <size_t N, typename ElemType, size_t M, size_t m>
RTree<N, ElemType, M, m>::~RTree() { root_pointer_.reset(); }

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::dimension() const {
  return N;
}

template <size_t N, typename ElemType, size_t M, size_t m>
size_t RTree<N, ElemType, M, m>::size() const {
  return entries;
}

template <size_t N, typename ElemType, size_t M, size_t m>
bool RTree<N, ElemType, M, m>::empty() const {
  return !entries;
}

template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::insert(const Rectangle<N> &box,
                                      const ElemType &value) {
  std::shared_ptr<Node> splitted_node = choose_leaf(root_pointer_, box, value);
  if (!splitted_node) {
    return;
  }
  // Split the root !
  std::shared_ptr<Node> new_root = std::make_shared<Node>();
  (*new_root)[0].child_pointer = root_pointer_;
  (new_root->size)++;
  adjust_tree(new_root, root_pointer_, splitted_node, &(*new_root)[0]);
  root_pointer_ = new_root;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_leaf(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      const ElemType &value) {
  if (!current_node->is_leaf()) {
    SpatialObject *entry;
    std::shared_ptr<Node> next_node = choose_node(current_node, box, entry);
    std::shared_ptr<Node> splitted_node = choose_leaf(next_node, box, value);
    return adjust_tree(current_node, next_node, splitted_node, entry);
  }
  SpatialObject new_entry;
  new_entry.box = box;
  new_entry.identifier = value;
  new_entry.child_pointer = nullptr;
  entries++;
  return current_node->insert(new_entry);
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::choose_node(const std::shared_ptr<Node> &current_node,
                                      const Rectangle<N> &box,
                                      SpatialObject *&entry) {
  Rectangle<N> enlarged_box = (*current_node)[0].box;
  enlarged_box.adjust(box);
  float minimum_area = (*current_node)[0].box.get_area();
  float minimum_enlargement = enlarged_box.get_area() - minimum_area;
  std::shared_ptr<Node> node = (*current_node)[0].child_pointer;
  float enlargement, area;

  entry = &(*current_node)[0];
  for (SpatialObject &current_entry : *current_node) {
    area = current_entry.box.get_area();
    enlarged_box = current_entry.box;
    enlarged_box.adjust(box);
    enlargement = enlarged_box.get_area() - area;

    if (enlargement < minimum_enlargement ||
        (enlargement == minimum_enlargement && area < minimum_area)) {
      minimum_enlargement = enlargement;
      minimum_area = area;
      node = current_entry.child_pointer;
      entry = &current_entry;
    }
  }
  return node;
}

template <size_t N, typename ElemType, size_t M, size_t m>
std::shared_ptr<typename RTree<N, ElemType, M, m>::Node>
RTree<N, ElemType, M, m>::adjust_tree(const std::shared_ptr<Node> &parent,
                                      const std::shared_ptr<Node> &first,
                                      const std::shared_ptr<Node> &second,
                                      SpatialObject *entry) {
  entry->box.reset();
  for (SpatialObject current_entry : *first) {
    entry->box.adjust(current_entry.box);
  }
  if (!second) {
    return nullptr;
  }
  SpatialObject new_entry;
  new_entry.box.reset();
  for (SpatialObject &current_entry : *second) {
    new_entry.box.adjust(current_entry.box);
  }
  new_entry.child_pointer = second;
  return parent->insert(new_entry);
}


template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::search(const std::shared_ptr<Node> &cur_node,
                                      const Rectangle<N> &box,
                                      std::vector<ElemType> &value){
  for (SpatialObject &obj : *cur_node) {
    if (overlaps(obj.box, box)) {
      if (cur_node->is_leaf()) {
        value.push_back(obj.identifier);
      }
      else {
        search(obj.child_pointer, box, value);
      }
    }
  }
}




template <size_t N, typename ElemType, size_t M, size_t m>
void RTree<N, ElemType, M, m>::print() {
  std::cout<<"Printing the tree ...";
  int lastLvl = -1;
  std::queue<std::pair<std::shared_ptr<Node>, int> > Q;
  Q.push(make_pair(root_pointer_, 0));
  while (!Q.empty()) {
    std::shared_ptr<Node> current = Q.front().first;
    int lvl = Q.front().second;
    Q.pop();
    if (lvl != lastLvl) {
      std::cout<<std::endl;
      lastLvl = lvl;
    }
    std::cout<<"{ | ";
    for (SpatialObject &obj : *current) {
      std::cout<<(obj.identifier.size() ? obj.identifier : "IN");
      std::cout<<":("<<obj.box[0].begin()<<","<<obj.box[0].end()<<") - ";
      std::cout<<"("<<obj.box[1].begin()<<","<<obj.box[1].end()<<") | ";
      if (!current->is_leaf()) {
        Q.push(make_pair(obj.child_pointer, lvl+1));
      }
    }
    std::cout<<"} ";    
  }
  std::cout<<"\n";
}

#endif  // SOURCE_RTREE_HPP_
