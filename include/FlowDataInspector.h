#pragma once

#include "FlowDataTypes.h"

#include <functional>

#include <QString>
#include <QStringList>

struct FlowDataScanResult {
    bool passed = false;
    QString errorMessage;
    QString warningMessage;
    QString sourceFolder;
    QStringList files;
    QStringList termCandidates;
    FlowDataEntry entry;
};

class FlowDataInspector {
public:
    using ProgressCallback = std::function<void(const QString& stage, int current, int total)>;
    using CancellationCallback = std::function<bool()>;

    // Accepts either a supported data directory or one supported data file.
    // This includes Tecplot ASCII (.dat) and binary (.plt) sources.
    static FlowDataScanResult scanFolder(const QString& folder,
                                         FlowDataKind dataKind,
                                         bool coordinateCheck,
                                         const ProgressCallback& progress = {},
                                         const CancellationCallback& shouldCancel = {});
    static FlowDataScanResult scanNek5000File(const QString& filePath,
                                              FlowDataKind dataKind,
                                              bool coordinateCheck,
                                              const ProgressCallback& progress = {},
                                              const CancellationCallback& shouldCancel = {});
    static FlowDataScanResult scanOpenFoamCase(const QString& casePath,
                                               const QString& caseType,
                                               FlowDataKind dataKind,
                                               bool coordinateCheck,
                                               bool startEndFieldsOnly,
                                               const ProgressCallback& progress = {},
                                               const CancellationCallback& shouldCancel = {});
    static FlowDataScanResult scanEnSightGoldCase(const QString& casePath,
                                                  FlowDataKind dataKind,
                                                  bool coordinateCheck,
                                                  bool commonFieldsOnly,
                                                  const ProgressCallback& progress = {},
                                                  const CancellationCallback& shouldCancel = {});
    static FlowDataScanResult scanEnSightGoldCases(const QStringList& casePaths,
                                                   FlowDataKind dataKind,
                                                   bool coordinateCheck,
                                                   bool commonFieldsOnly,
                                                   const ProgressCallback& progress = {},
                                                   const CancellationCallback& shouldCancel = {});
    static bool convertOpenFoamCaseToPvd(const QString& casePath,
                                         const QString& caseType,
                                         const QString& destinationFolder,
                                         const QString& pvdBaseName,
                                         bool startEndFieldsOnly,
                                         const QStringList& allowedFields,
                                         QString* errorMessage = nullptr,
                                         const ProgressCallback& progress = {},
                                         const CancellationCallback& shouldCancel = {});
    static bool convertEnSightGoldCaseToPvd(const QString& casePath,
                                            const QString& destinationFolder,
                                            const QString& pvdBaseName,
                                            bool commonFieldsOnly,
                                            const QStringList& allowedFields,
                                            QString* errorMessage = nullptr,
                                            const ProgressCallback& progress = {},
                                            const CancellationCallback& shouldCancel = {});
    static bool convertEnSightGoldCasesToPvd(const QStringList& casePaths,
                                             const QString& destinationFolder,
                                             const QString& pvdBaseName,
                                             bool commonFieldsOnly,
                                             const QStringList& allowedFields,
                                             QString* errorMessage = nullptr,
                                             const ProgressCallback& progress = {},
                                             const CancellationCallback& shouldCancel = {});
};
