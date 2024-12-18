//
// Created by Grzegorz Bokota on 24.10.24.
//

#ifndef PARTSEGCORE_TRIANGULATION_DEBUG_UTIL_
#define PARTSEGCORE_TRIANGULATION_DEBUG_UTIL_

#include <iostream>

namespace partsegcore {

template <typename T>
void print_set(std::ostream &o, const T &s, const std::string &end = "\n") {
  bool first = true;
  o << "{";
  for (const auto &el : s) {
    if (!first) {
      o << ", ";
    }
    o << el;
    first = false;
  }
  o << "}" << end;
}
template <typename T>
void print_vector(std::ostream &o, const T &s, const std::string &end = "\n") {
  bool first = true;
  o << "[";
  for (const auto &el : s) {
    if (!first) {
      o << ", ";
    }
    o << el;
    first = false;
  }
  o << "]" << end;
}
template <typename T>
void print_map(std::ostream &o, const T &s, const std::string &end = "\n") {
  bool first = true;
  o << "{";
  for (const auto &el : s) {
    if (!first) {
      o << ", ";
    }
    o << el.first << ": ";
    first = false;
    print_set(o, el.second, "");
  }
  o << "}" << end;
}

template <typename T1, typename T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) {
  os << "(" << p.first << ", " << p.second << ")";
  return os;
}
}  // namespace partsegcore

#endif  // PARTSEGCORE_TRIANGULATION_DEBUG_UTIL_
