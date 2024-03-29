#ifndef HELPER_F_H
#define HELPER_F_H
#include <algorithm>
#include <numeric>
#include <limits>
#include <list>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>
#include "RcppArmadillo.h"
using namespace Rcpp;

//replicate elements in a vector
template<typename T>
std::vector<T> replicate(const T& val, size_t n) {
    return std::vector<T>(n, val);
}
std::vector<size_t> gen(const std::vector<std::vector<double>>& Lu);
//equivalent to R's unlist() function
template<typename T>
std::vector<T> unlist(const std::vector<std::vector<T>>& vec)
{
  int size = 0;
  for(auto it = vec.begin(); it != vec.end();++it)
  {
    size += it->size();
  }
  std::vector<T> unrolled;
  unrolled.reserve(size);
  for (const auto& inner_vec : vec) {
      for (const auto& item : inner_vec) {
          unrolled.push_back(item);
      }
  }
  return unrolled;
}

//moving average filter
void filter(std::vector<double>& values, const double& a);

std::vector<double> make_bins(std::set<double>&, int, bool);

//cut(ydata$.index, breaks = bins, include.lowest = TRUE)
std::vector<std::pair<double,double>> cut(const std::vector<double>& y_data_index, const std::set<double> breaks);


std::vector<size_t> generateIdVec(const std::vector<std::vector<double>>& Ly);
#endif