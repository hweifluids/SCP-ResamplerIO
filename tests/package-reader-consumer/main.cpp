#include <streamcenterplus/io/DataSetReader.h>

int main() {
    using streamcenterplus::io::DataSetFrameKind;
    return streamcenterplus::io::ToString(DataSetFrameKind::Nek5000TimeStep) ==
                   "nek5000_time_step"
               ? 0
               : 1;
}
