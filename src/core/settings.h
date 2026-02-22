#pragma once
#include <string>

#if __cplusplus > 201703L
constexpr std::string OUTPUT_DIR = "../output";
#else
const std::string OUTPUT_DIR = "../output";
#endif
