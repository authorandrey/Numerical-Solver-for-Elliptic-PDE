#pragma once
#include <string_view>

#if __cplusplus > 201703L
constexpr std::string OUTPUT_DIR = "../output";
#else
const std::string OUTPUT_DIR = "../output";
#endif