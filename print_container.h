#ifndef PRINT_CONTAINER_H
#define PRINT_CONTAINER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <map>

template <typename T, typename S>
std::ostream& operator<<(std::ostream &os, std::pair<T,S> p)
{
  os << p.first << "=>" << p.second;
  return os;
}


template <typename T>
void print_container(std::ostream &os, T it, T it_end,
		     const char *pre, const char *delim, const char *post)
{
  if (it == it_end)
    return;

  os << pre;
  os << *it;
  ++it;
  while (it != it_end)
  {
    os << delim << *it;
    ++it;
  }
  os << post;
}

template<typename T>
void print_container(std::ostream &os, T c,
		     const char *pre, const char *delim, const char *post)
{
  print_container(os, c.begin(), c.end(), pre, delim, post);
}

#endif
