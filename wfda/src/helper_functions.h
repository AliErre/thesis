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

template<typename T>
std::vector<int> generateIdVec(const std::vector<std::vector<T>>& Ly) {
    std::vector<int> id_vec;
    for (size_t i = 0; i < Ly.size(); ++i) {
        id_vec.reserve(id_vec.size() + Ly[i].size());
        for (size_t j = 0; j < Ly[i].size(); ++j) {
            id_vec.push_back(static_cast<int>(i)); //(i+1)? to id_vec
        }
    }
    id_vec.shrink_to_fit();//non necessario credo
    return id_vec;
}

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


#endif