#ifndef RECONSTRUCTIONFACTORY_H
#define RECONSTRUCTIONFACTORY_H

#include "reconstruction.h" // contains the abstract base class and its subclasses
#include <memory>
#include <string>

std::unique_ptr<ReconstructionBase> reconstructionFactory(const std::string&, const NumericMatrix&, 
                                                          double, int, NumericVector&,
                                                          int, int);

#endif // RECONSTRUCTIONFACTORY_H