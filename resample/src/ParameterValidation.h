#pragma once

#include "streamcenterplus_resample/Core.h"

namespace streamcenterplus::resample {

// Normalizes defaults and rejects values that cannot be executed safely.
// Kept internal so both deck loading and direct RunSolver callers share the
// same validation boundary before any output is created or removed.
void ValidateParameters(Parameters* params);

}  // namespace streamcenterplus::resample
