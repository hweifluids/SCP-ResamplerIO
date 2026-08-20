#pragma once

#include "streamcenterplus_resample/Core.h"

namespace streamcenterplus::resample {

// Kriging solves a dense (neighbors + 1)-square system for every output point.
// Keep the public deck boundary well below sizes that could overflow indexing or
// turn an accidental input value into an impractical allocation.
inline constexpr int kMaximumKrigingNeighbors = 1024;

// Normalizes defaults and rejects values that cannot be executed safely.
// Kept internal so both deck loading and direct RunSolver callers share the
// same validation boundary before any output is created or removed.
void ValidateParameters(Parameters* params);

}  // namespace streamcenterplus::resample
