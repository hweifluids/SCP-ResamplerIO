#include "streamcenterplus_resample/Core.h"

#include <streamcenterplus/SolverInfo.h>

int main(int argc, char** argv) {
    if (streamcenterplus::HandleSolverInfoRequest(argc, argv)) {
        return 0;
    }
    return streamcenterplus::resample::RunCli(argc, argv);
}
