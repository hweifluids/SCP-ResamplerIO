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
    QVector<FlowDataSourceFingerprint> sourceInventory;
    FlowDataEntry entry;
};

class FlowDataInspector {
public:
    using ProgressCallback = std::function<void(const QString& stage, int current, int total)>;
    using CancellationCallback = std::function<bool()>;

    // Accepts either a supported data directory or one supported data file.
    // This includes Tecplot ASCII (.dat) and binary (.plt) sources. Resampled
    // sources always validate every frame's uniform grid, dimensions, and
    // numeric point/cell array schema; coordinateCheck additionally compares
    // the frame coordinates.
    static FlowDataScanResult scanFolder(const QString& folder,
                                         FlowDataKind dataKind,
                                         bool coordinateCheck,
                                         const ProgressCallback& progress = {},
                                         const CancellationCallback& shouldCancel = {});
    // Validates converter staging output and derives dt from its authoritative
    // PVD time axis in addition to the ordinary VTK metadata scan.
    static FlowDataScanResult scanConvertedPvdFolder(const QString& folder,
                                                     FlowDataKind dataKind,
                                                     bool coordinateCheck,
                                                     const ProgressCallback& progress = {},
                                                     const CancellationCallback& shouldCancel = {});
    // Reads Native Nek5000 through the shared canonical element-local GLL
    // reader.  It deliberately does not construct legacy linear GLL subcells
    // while probing frames, coordinates, or fields.
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
    // Accepts .case and Fluent .encas descriptors. A one-step geometry time
    // set is reusable with transient variables. Decimal export precision may
    // support an explicitly warned effective dt; original native frame times
    // remain unchanged through snapshots and PVD conversion.
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
    static bool buildEnSightGoldSourceInventory(
        const QStringList& casePaths,
        QVector<FlowDataSourceFingerprint>* inventory,
        QString* errorMessage = nullptr,
        const ProgressCallback& progress = {},
        const CancellationCallback& shouldCancel = {});
    // Copies the exact validated dependency inventory into a private caller-
    // owned directory, then emits normalized case files that reference only
    // immutable snapshot copies.  The converter must use the returned paths,
    // never the original case paths.
    static bool materializeEnSightGoldSourceSnapshot(
        const QStringList& casePaths,
        const QVector<FlowDataSourceFingerprint>& expectedInventory,
        const QString& snapshotRoot,
        QStringList* snapshotCasePaths,
        QString* errorMessage = nullptr,
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
