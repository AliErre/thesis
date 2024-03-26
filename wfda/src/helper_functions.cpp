#include "helper_functions.h"

void filter(std::vector<double>& values, const double& a) {
    for (size_t i = 0; i < values.size() - 1; ++i) { // Skip the last element
        values[i] = a * values[i] + (1 - a) * values[i + 1];
    }
    values.pop_back();
}

std::vector<double> make_bins(std::set<double>& bins, int maxbins, bool f) {
  
  double minValue = *bins.begin(); //sono ordinati in un set
  double maxValue = *bins.rbegin();
  
  if(f){
    std::vector<double> binvalues;
    binvalues.reserve(maxbins + 1);
    minValue *= (1 - 0.001 * ((minValue > 0) - (minValue < 0)));//(min > 0) - (min < 0) restituisce il segno
    maxValue *= (1 + 0.001 * ((maxValue > 0) - (maxValue < 0)));

    //sequence
    double step = (maxValue - minValue) / maxbins;
    for (int i = 0; i <= maxbins; ++i) {
        binvalues.push_back(minValue + i * step);
    }
    bins.clear();
    bins.insert(binvalues.begin(),binvalues.end());

    filter(binvalues, 0.5);//memory managed automatically
    return binvalues;
  }else{
    std::vector<double> binvalues(bins.begin(), bins.end());
    
    double start = (1 - 0.001 * (minValue > 0.0 ? 1 : (minValue < 0.0 ? -1 : 0.0))) * minValue;
    double end = (1 + 0.001 * (maxValue > 0.0 ? 1 : (maxValue < 0.0 ? -1 : 0.0))) * maxValue;
    if(start == 0) start = -0.001;
    if (end == 0) end = 0.001;

    bins.erase(std::prev(bins.end()));//cant do iteartor - 1 for sets cause theyre not contiguous elements
    std::set<double> temp = bins;
    bins.clear();
    bins.insert(start);
    bins.insert(temp.begin(),temp.end());
    bins.insert(end);

    return binvalues;
  }

}

//cut(ydata$.index, breaks = bins, include.lowest = TRUE)
std::vector<std::pair<double,double>> cut(const std::vector<double>& y_data_index, const std::set<double> breaks)
{
  std::vector<std::pair<double,double>> pairs;
  pairs.reserve(y_data_index.size());//avrò tante classi quanti ydata.size()

  for(auto& elem:y_data_index)
  {
    auto it = breaks.lower_bound(elem);
    if(it == breaks.begin())
    {
      auto next = std::next(it);
      pairs.emplace_back(*it,*next);
    }
    else{
      auto prev = std::prev(it);
      pairs.emplace_back(*prev,*it);//dovrebbe andare bene anche per il caso it == std::prev(breaks.end());  
    }
  }

  pairs.shrink_to_fit();
  return pairs;
}
