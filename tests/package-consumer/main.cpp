#include <FlowDataInspector.h>
#include <FlowDataTypes.h>
#include <streamcenterplus_resample/Core.h>

int main() {
    FlowDataEntry entry;
    entry.dataKind = FlowDataKind::Raw;
    entry.source = FlowDataSource::Link;
    streamcenterplus::resample::FrameInfo frame;
    frame.kind = streamcenterplus::resample::DataSetFrameKind::Nek5000TimeStep;
    return streamcenterplus::resample::ToString(frame.kind) == "nek5000_time_step" ? 0 : 1;
}
