#ifndef RECONSTRUCTIONFACTORY_H
#define RECONSTRUCTIONFACTORY_H

#include "reconstruction.h"
#include <memory>
#include <string>

//std::unique_ptr<ReconstructionBase> reconstructionFactory(const std::string&, const NumericMatrix&);
ReconstructionBase* reconstructionFactory(const std::string& id, const NumericMatrix& Y);
#endif // RECONSTRUCTIONFACTORY_H