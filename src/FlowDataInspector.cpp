#include "FlowDataInspector.h"
#include "FlowDataProtocolLimits.h"
#include "streamcenterplus_resample/Core.h"

#include <streamcenterplus/io/DataSetReader.h>
#include <streamcenterplus/VtkHdfIO.h>

#include <vtkAbstractArray.h>
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkGenericEnSightReader.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkOpenFOAMReader.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPOpenFOAMReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <QDate>
#include <QCryptographicHash>
#include <QByteArray>
#include <QDir>
#include <QDirIterator>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSet>
#include <QStorageInfo>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <limits>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#ifdef Q_OS_WIN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {

bool isFileSystemReparsePoint(const QString& path) {
#ifdef Q_OS_WIN
    const std::wstring nativePath = QDir::toNativeSeparators(path).toStdWString();
    const DWORD attributes = GetFileAttributesW(nativePath.c_str());
    return attributes != INVALID_FILE_ATTRIBUTES
        && (attributes & FILE_ATTRIBUTE_REPARSE_POINT) != 0;
#else
    return QFileInfo(path).isSymLink();
#endif
}

bool directoryChainContainsReparsePoint(const QString& directoryPath) {
    QString current = QFileInfo(directoryPath).absoluteFilePath();
    while (!current.isEmpty()) {
        const QFileInfo info(current);
        if (!info.exists() || !info.isDir() || info.isSymLink()
            || isFileSystemReparsePoint(current)) {
            return true;
        }
        const QString parent = info.absoluteDir().absolutePath();
        if (parent == current) {
            break;
        }
        current = parent;
    }
    return false;
}

bool makeSnapshotFileOwnerReadOnly(const QString& path, QString* errorMessage) {
#ifdef Q_OS_WIN
    const std::wstring nativePath = QDir::toNativeSeparators(path).toStdWString();
    const DWORD attributes = GetFileAttributesW(nativePath.c_str());
    if (attributes == INVALID_FILE_ATTRIBUTES
        || (attributes & (FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_REPARSE_POINT)) != 0
        || !SetFileAttributesW(nativePath.c_str(), attributes | FILE_ATTRIBUTE_READONLY)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot make EnSight snapshot file read-only: %1")
                                .arg(path);
        }
        return false;
    }
    const DWORD verifiedAttributes = GetFileAttributesW(nativePath.c_str());
    if (verifiedAttributes == INVALID_FILE_ATTRIBUTES
        || (verifiedAttributes & FILE_ATTRIBUTE_READONLY) == 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight snapshot file did not remain read-only: %1")
                                .arg(path);
        }
        return false;
    }
    return true;
#else
    if (!QFile::setPermissions(path, QFileDevice::ReadOwner)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot make EnSight snapshot file owner-read-only: %1")
                                .arg(path);
        }
        return false;
    }
    return true;
#endif
}

QString normalizedPath(const QString& path) {
    return QDir::fromNativeSeparators(path);
}

std::filesystem::path fsPathFromQString(const QString& path) {
#ifdef Q_OS_WIN
    return std::filesystem::path(QDir::toNativeSeparators(path).toStdWString());
#else
    const QByteArray encoded = QFile::encodeName(QDir::toNativeSeparators(path));
    return std::filesystem::path(encoded.constData());
#endif
}

QString qStringFromFsPath(const std::filesystem::path& path) {
#ifdef Q_OS_WIN
    return normalizedPath(QString::fromStdWString(path.wstring()));
#else
    return normalizedPath(QFile::decodeName(path.native().c_str()));
#endif
}

QString realVtkPath(const QString& path) {
    try {
        return qStringFromFsPath(
            streamcenterplus::vtkhdf::RealPath(fsPathFromQString(path)));
    } catch (const std::filesystem::filesystem_error&) {
        return {};
    }
}

QString vtkExtensionSuffix(const QString& path) {
    std::string ext;
    try {
        ext = streamcenterplus::vtkhdf::Extension(fsPathFromQString(path));
    } catch (const std::filesystem::filesystem_error&) {
        return {};
    }
    if (!ext.empty() && ext.front() == '.') {
        ext.erase(ext.begin());
    }
    return QString::fromStdString(ext);
}

QString formatDouble(double value) {
    return QString::number(value, 'g', 15);
}

QString formatMachineDouble(double value) {
    return QString::number(value, 'g', std::numeric_limits<double>::max_digits10);
}

QString formatRange(double minValue, double maxValue) {
    return QStringLiteral("[%1, %2]").arg(formatDouble(minValue), formatDouble(maxValue));
}

bool supportedVtkExtension(const QString& suffix, bool includeNek5000) {
    const QString ext = suffix.toLower();
    return ext == QStringLiteral("vtk")
        || ext == QStringLiteral("vti")
        || ext == QStringLiteral("vtr")
        || ext == QStringLiteral("vts")
        || ext == QStringLiteral("vtu")
        || ext == QStringLiteral("vtm")
        || ext == QStringLiteral("vtkhdf")
        || ext == QStringLiteral("dat")
        || ext == QStringLiteral("plt")
        || (includeNek5000 &&
            (ext == QStringLiteral("nek5000") || ext == QStringLiteral("nek")));
}

bool isTecplotExtension(const QString& suffix) {
    const QString ext = suffix.toLower();
    return ext == QStringLiteral("dat") || ext == QStringLiteral("plt");
}

bool isNek5000DescriptorExtension(const QString& suffix) {
    const QString ext = suffix.toLower();
    return ext == QStringLiteral("nek5000") || ext == QStringLiteral("nek");
}

QString normalizedOpenFoamCaseType(QString caseType) {
    caseType = caseType.trimmed().toLower();
    if (caseType == QStringLiteral("reconstructed")
        || caseType.contains(QStringLiteral("reconstruct"))) {
        return QStringLiteral("reconstructed");
    }
    if (caseType == QStringLiteral("decomposed")
        || caseType.contains(QStringLiteral("decompos"))) {
        return QStringLiteral("decomposed");
    }
    return {};
}

QString openFoamFileType(const QString& caseType) {
    const QString normalized = normalizedOpenFoamCaseType(caseType);
    return normalized == QStringLiteral("decomposed")
        ? QStringLiteral("openfoam decomposed case (pvd+vtu)")
        : QStringLiteral("openfoam reconstructed case (pvd+vtu)");
}

QString openFoamCaseRoot(const QString& casePath) {
    return normalizedPath(QFileInfo(casePath).absolutePath());
}

QString enSightGoldFileType() {
    return QStringLiteral("ensight gold case (pvd+vtu)");
}

QString stripOpenFoamComments(QString text) {
    text.remove(QRegularExpression(QStringLiteral(R"(/\*.*?\*/)"),
                                   QRegularExpression::DotMatchesEverythingOption));
    text.remove(QRegularExpression(QStringLiteral(R"(//[^\r\n]*)")));
    return text;
}

bool readOpenFoamWriteInterval(const QString& casePath, double* writeInterval, QString* errorMessage) {
    if (writeInterval != nullptr) {
        *writeInterval = 0.0;
    }
    const QString controlDictPath = normalizedPath(QDir(openFoamCaseRoot(casePath)).absoluteFilePath(
        QStringLiteral("system/controlDict")));
    QFile file(controlDictPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot read OpenFOAM system/controlDict: %1").arg(controlDictPath);
        }
        return false;
    }

    const QString text = stripOpenFoamComments(QString::fromUtf8(file.readAll()));
    const QRegularExpression writeIntervalRegex(QStringLiteral(R"(\bwriteInterval\s+([^;\s]+)\s*;)"));
    const QRegularExpressionMatch match = writeIntervalRegex.match(text);
    if (!match.hasMatch()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("OpenFOAM system/controlDict does not define writeInterval: %1").arg(controlDictPath);
        }
        return false;
    }

    bool ok = false;
    const double value = match.captured(1).trimmed().toDouble(&ok);
    if (!ok || !std::isfinite(value) || value <= 0.0) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("OpenFOAM writeInterval must be a positive numeric value: %1").arg(match.captured(1).trimmed());
        }
        return false;
    }
    if (writeInterval != nullptr) {
        *writeInterval = value;
    }
    return true;
}

qlonglong lastNumericTag(const QString& fileName) {
    static const QRegularExpression numberRegex(QStringLiteral("(\\d+)"));
    qlonglong tag = -1;
    const auto matches = numberRegex.globalMatch(fileName);
    auto it = matches;
    while (it.hasNext()) {
        const QRegularExpressionMatch match = it.next();
        bool ok = false;
        const qlonglong value = match.captured(1).toLongLong(&ok);
        if (ok) {
            tag = value;
        }
    }
    return tag;
}

QStringList collectVtkFiles(const QString& sourcePath, bool includeNek5000) {
    const QFileInfo sourceInfo(sourcePath);
    QFileInfoList entries;
    if (sourceInfo.exists() && sourceInfo.isFile() && sourceInfo.isReadable()) {
        entries.append(sourceInfo);
    } else {
        QDir dir(sourcePath);
        entries = dir.entryInfoList(QDir::Files | QDir::NoDotAndDotDot | QDir::Readable,
                                    QDir::Name | QDir::IgnoreCase);
    }
    QList<QPair<QPair<qlonglong, QString>, QString>> tagged;
    for (const QFileInfo& info : entries) {
        if (!supportedVtkExtension(info.suffix(), includeNek5000)) {
            continue;
        }
        const std::filesystem::path filePath = fsPathFromQString(info.absoluteFilePath());
        std::vector<std::filesystem::path> framePaths;
        try {
            framePaths = streamcenterplus::vtkhdf::ExpandTimeSteps(filePath);
        } catch (const std::exception&) {
            framePaths = {filePath};
        }
        for (const std::filesystem::path& framePath : framePaths) {
            const QString frameText = qStringFromFsPath(framePath);
            tagged.push_back(qMakePair(qMakePair(lastNumericTag(info.fileName()), frameText), frameText));
        }
    }
    std::sort(tagged.begin(), tagged.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first.first != rhs.first.first) {
            return lhs.first.first < rhs.first.first;
        }
        return lhs.first.second.compare(rhs.first.second, Qt::CaseInsensitive) < 0;
    });

    QStringList files;
    files.reserve(tagged.size());
    for (const auto& entry : tagged) {
        files.append(entry.second);
    }
    return files;
}

vtkSmartPointer<vtkOpenFOAMReader> createOpenFoamReader(const QString& caseType) {
    auto reader = vtkSmartPointer<vtkPOpenFOAMReader>::New();
    reader->SetCaseType(normalizedOpenFoamCaseType(caseType) == QStringLiteral("decomposed")
                            ? vtkPOpenFOAMReader::DECOMPOSED_CASE
                            : vtkPOpenFOAMReader::RECONSTRUCTED_CASE);
    return reader;
}

bool prepareOpenFoamReader(vtkOpenFOAMReader* reader,
                           const QString& casePath,
                           const QString& caseType,
                           QString* errorMessage) {
    if (reader == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: OpenFOAM reader was not created.");
        }
        return false;
    }

    const QString normalizedCaseType = normalizedOpenFoamCaseType(caseType);
    if (normalizedCaseType.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Case type must be Reconstructed or Decomposed.");
        }
        return false;
    }

    const QFileInfo caseInfo(casePath);
    if (!caseInfo.exists() || !caseInfo.isFile()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("OpenFOAM foam file does not exist: %1").arg(casePath);
        }
        return false;
    }
    if (caseInfo.suffix().compare(QStringLiteral("foam"), Qt::CaseInsensitive) != 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("OpenFOAM import expects a *.foam file.");
        }
        return false;
    }

    const QByteArray fileName = QFile::encodeName(caseInfo.absoluteFilePath());
    reader->SetFileName(fileName.constData());
    reader->SetCreateCellToPoint(false);
    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->DisableAllPatchArrays();
    reader->EnableAllCellArrays();
    reader->EnableAllPointArrays();
    reader->UpdateInformation();

    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->EnableAllCellArrays();
    reader->EnableAllPointArrays();
    reader->DisableAllPatchArrays();

    bool enabledInternalMesh = false;
    for (int index = 0; index < reader->GetNumberOfPatchArrays(); ++index) {
        const char* rawName = reader->GetPatchArrayName(index);
        const QString name = QString::fromUtf8(rawName == nullptr ? "" : rawName).trimmed();
        if (name.compare(QStringLiteral("internalMesh"), Qt::CaseInsensitive) == 0
            || name.contains(QStringLiteral("internal"), Qt::CaseInsensitive)) {
            reader->SetPatchArrayStatus(rawName, 1);
            enabledInternalMesh = true;
        }
    }
    Q_UNUSED(enabledInternalMesh);
    return true;
}

QVector<double> openFoamTimeValues(vtkOpenFOAMReader* reader) {
    QVector<double> values;
    if (reader == nullptr) {
        return values;
    }
    vtkDoubleArray* timeValues = reader->GetTimeValues();
    if (timeValues == nullptr) {
        return values;
    }
    values.reserve(static_cast<int>(std::min<vtkIdType>(timeValues->GetNumberOfTuples(),
                                                        std::numeric_limits<int>::max())));
    for (vtkIdType index = 0; index < timeValues->GetNumberOfTuples(); ++index) {
        const double value = timeValues->GetValue(index);
        if (std::isfinite(value)) {
            values.push_back(value);
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

void collectVolumeDataSets(vtkDataObject* object, QList<vtkSmartPointer<vtkDataSet>>* dataSets) {
    if (object == nullptr || dataSets == nullptr) {
        return;
    }
    if (auto* dataSet = vtkDataSet::SafeDownCast(object)) {
        if (vtkPolyData::SafeDownCast(dataSet) == nullptr
            && dataSet->GetNumberOfPoints() > 0
            && dataSet->GetNumberOfCells() > 0) {
            dataSets->append(dataSet);
        }
        return;
    }
    if (auto* composite = vtkCompositeDataSet::SafeDownCast(object)) {
        vtkSmartPointer<vtkCompositeDataIterator> iterator;
        iterator.TakeReference(composite->NewIterator());
        if (iterator == nullptr) {
            return;
        }
        iterator->SkipEmptyNodesOn();
        for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal(); iterator->GoToNextItem()) {
            collectVolumeDataSets(iterator->GetCurrentDataObject(), dataSets);
        }
    }
}

vtkSmartPointer<vtkDataSet> extractOpenFoamVolume(vtkDataObject* output, QString* errorMessage) {
    QList<vtkSmartPointer<vtkDataSet>> dataSets;
    collectVolumeDataSets(output, &dataSets);
    if (dataSets.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The OpenFOAM reader did not produce any volume dataset.");
        }
        return nullptr;
    }

    auto append = vtkSmartPointer<vtkAppendFilter>::New();
    append->MergePointsOn();
    for (const vtkSmartPointer<vtkDataSet>& dataSet : std::as_const(dataSets)) {
        append->AddInputData(dataSet);
    }
    append->Update();

    auto retained = vtkSmartPointer<vtkUnstructuredGrid>::New();
    retained->ShallowCopy(append->GetOutput());
    if (retained->GetNumberOfPoints() <= 0 || retained->GetNumberOfCells() <= 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The OpenFOAM volume dataset is empty.");
        }
        return nullptr;
    }
    return retained;
}

vtkSmartPointer<vtkDataSet> readOpenFoamVolumeAtTime(vtkOpenFOAMReader* reader,
                                                     double timeValue,
                                                     QString* errorMessage) {
    if (reader == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: OpenFOAM reader is null.");
        }
        return nullptr;
    }
    reader->SetTimeValue(timeValue);
    reader->UpdateTimeStep(timeValue);
    vtkDataObject* output = reader->GetOutputDataObject(0);
    return extractOpenFoamVolume(output, errorMessage);
}

bool prepareEnSightGoldReader(vtkGenericEnSightReader* reader,
                              const QString& casePath,
                              QString* errorMessage) {
    if (reader == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight Gold reader was not created.");
        }
        return false;
    }

    const QFileInfo caseInfo(casePath);
    if (!caseInfo.exists() || !caseInfo.isFile() || caseInfo.isSymLink()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("EnSight Gold case file does not exist or is not a regular non-link file: %1")
                                .arg(casePath);
        }
        return false;
    }
    if (caseInfo.suffix().compare(QStringLiteral("case"), Qt::CaseInsensitive) != 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold import expects a *.case file.");
        }
        return false;
    }

    const QByteArray fileName = QFile::encodeName(caseInfo.absoluteFilePath());
    reader->SetCaseFileName(fileName.constData());
    reader->ReadAllVariablesOn();
    reader->UpdateInformation();

    const int version = reader->GetEnSightVersion();
    if (version != vtkGenericEnSightReader::ENSIGHT_GOLD
        && version != vtkGenericEnSightReader::ENSIGHT_GOLD_BINARY) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The selected *.case file is not an EnSight Gold case.");
        }
        return false;
    }

    reader->ReadAllVariablesOn();
    for (int index = 0; index < reader->GetNumberOfPointArrays(); ++index) {
        const char* name = reader->GetPointArrayName(index);
        if (name != nullptr) {
            reader->SetPointArrayStatus(name, 1);
        }
    }
    for (int index = 0; index < reader->GetNumberOfCellArrays(); ++index) {
        const char* name = reader->GetCellArrayName(index);
        if (name != nullptr) {
            reader->SetCellArrayStatus(name, 1);
        }
    }
    return true;
}

QVector<double> enSightGoldTimeValues(vtkGenericEnSightReader* reader) {
    QVector<double> values;
    if (reader == nullptr) {
        return values;
    }

    vtkInformation* information = reader->GetOutputInformation(0);
    if (information != nullptr) {
        const int count = information->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        double* rawValues = information->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        for (int index = 0; index < count; ++index) {
            const double value = rawValues[index];
            if (std::isfinite(value)) {
                values.push_back(value);
            }
        }
    }

    if (values.isEmpty()) {
        vtkDataArrayCollection* timeSets = reader->GetTimeSets();
        if (timeSets != nullptr) {
            for (int setIndex = 0; setIndex < timeSets->GetNumberOfItems(); ++setIndex) {
                vtkDataArray* array = timeSets->GetItem(setIndex);
                if (array == nullptr) {
                    continue;
                }
                for (vtkIdType tuple = 0; tuple < array->GetNumberOfTuples(); ++tuple) {
                    const double value = array->GetComponent(tuple, 0);
                    if (std::isfinite(value)) {
                        values.push_back(value);
                    }
                }
            }
        }
    }

    if (values.isEmpty()) {
        values.push_back(0.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

double inferredTimeInterval(const QVector<double>& timeValues) {
    QVector<double> positiveSteps;
    for (int index = 1; index < timeValues.size(); ++index) {
        const double step = timeValues.at(index) - timeValues.at(index - 1);
        if (std::isfinite(step) && step > 0.0) {
            positiveSteps.push_back(step);
        }
    }
    if (positiveSteps.isEmpty()) {
        return 1.0;
    }
    std::sort(positiveSteps.begin(), positiveSteps.end());
    return positiveSteps.at(positiveSteps.size() / 2);
}

bool uniformTimeInterval(const QVector<double>& timeValues,
                         double* interval,
                         QString* errorMessage) {
    if (interval == nullptr) {
        return false;
    }
    if (timeValues.size() <= 1) {
        *interval = 1.0;
        return true;
    }
    const double firstStep = timeValues.at(1) - timeValues.at(0);
    if (!std::isfinite(firstStep) || firstStep <= 0.0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight timesteps must be finite and strictly increasing.");
        }
        return false;
    }
    for (qsizetype index = 2; index < timeValues.size(); ++index) {
        const double current = timeValues.at(index);
        const double previous = timeValues.at(index - 1);
        const double expected = std::fma(static_cast<double>(index),
                                         firstStep, timeValues.first());
        const double magnitude = std::max(std::abs(current), std::abs(expected));
        const double tolerance = std::abs(firstStep)
                * static_cast<double>(index) * 1.0e-10
            + std::numeric_limits<double>::epsilon() * magnitude * 8.0
            + std::numeric_limits<double>::denorm_min() * 16.0;
        if (!std::isfinite(current) || !(current > previous)
            || !std::isfinite(expected)
            || std::abs(current - expected) > tolerance) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight uses a non-uniform time axis that cannot be represented by one dt value.");
            }
            return false;
        }
    }
    *interval = firstStep;
    return true;
}

vtkSmartPointer<vtkDataSet> extractEnSightGoldVolume(vtkDataObject* output, QString* errorMessage) {
    QList<vtkSmartPointer<vtkDataSet>> dataSets;
    collectVolumeDataSets(output, &dataSets);
    if (dataSets.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The EnSight Gold reader did not produce any volume dataset.");
        }
        return nullptr;
    }

    auto append = vtkSmartPointer<vtkAppendFilter>::New();
    append->MergePointsOn();
    for (const vtkSmartPointer<vtkDataSet>& dataSet : std::as_const(dataSets)) {
        append->AddInputData(dataSet);
    }
    append->Update();

    auto retained = vtkSmartPointer<vtkUnstructuredGrid>::New();
    retained->ShallowCopy(append->GetOutput());
    if (retained->GetNumberOfPoints() <= 0 || retained->GetNumberOfCells() <= 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The EnSight Gold volume dataset is empty.");
        }
        return nullptr;
    }
    return retained;
}

vtkSmartPointer<vtkDataSet> readEnSightGoldVolumeAtTime(vtkGenericEnSightReader* reader,
                                                        double timeValue,
                                                        QString* errorMessage) {
    if (reader == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight Gold reader is null.");
        }
        return nullptr;
    }
    reader->SetTimeValue(timeValue);
    reader->UpdateTimeStep(timeValue);
    vtkDataObject* output = reader->GetOutputDataObject(0);
    return extractEnSightGoldVolume(output, errorMessage);
}

vtkSmartPointer<vtkDataSet> readEnSightGoldVolumeFromCaseAtTime(const QString& casePath,
                                                                double timeValue,
                                                                QString* errorMessage) {
    auto reader = vtkSmartPointer<vtkGenericEnSightReader>::New();
    if (!prepareEnSightGoldReader(reader, casePath, errorMessage)) {
        return nullptr;
    }
    return readEnSightGoldVolumeAtTime(reader, timeValue, errorMessage);
}

bool writeUnstructuredFrame(vtkDataSet* dataSet, const QString& filePath, QString* errorMessage) {
    if (dataSet == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: no VTK dataset to write.");
        }
        return false;
    }

    auto append = vtkSmartPointer<vtkAppendFilter>::New();
    append->MergePointsOn();
    append->AddInputData(dataSet);
    append->Update();

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(QFile::encodeName(filePath).constData());
    writer->SetInputData(append->GetOutput());
    writer->SetDataModeToBinary();
    if (writer->Write() != 1) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Failed to write VTK frame: %1").arg(filePath);
        }
        return false;
    }
    return true;
}

bool writePvdFile(const QString& pvdPath,
                  const QVector<double>& timeValues,
                  const QStringList& relativeFramePaths,
                  QString* errorMessage) {
    if (timeValues.size() != relativeFramePaths.size() || timeValues.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: inconsistent PVD frame list.");
        }
        return false;
    }

    QSaveFile file(pvdPath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot write PVD file: %1").arg(file.errorString());
        }
        return false;
    }
    QXmlStreamWriter xml(&file);
    xml.setAutoFormatting(true);
    xml.writeStartDocument();
    xml.writeStartElement(QStringLiteral("VTKFile"));
    xml.writeAttribute(QStringLiteral("type"), QStringLiteral("Collection"));
    xml.writeAttribute(QStringLiteral("version"), QStringLiteral("0.1"));
    xml.writeAttribute(QStringLiteral("byte_order"), QStringLiteral("LittleEndian"));
    xml.writeStartElement(QStringLiteral("Collection"));
    for (int index = 0; index < timeValues.size(); ++index) {
        xml.writeEmptyElement(QStringLiteral("DataSet"));
        xml.writeAttribute(QStringLiteral("timestep"),
                           formatMachineDouble(timeValues.at(index)));
        xml.writeAttribute(QStringLiteral("group"), QString());
        xml.writeAttribute(QStringLiteral("part"), QStringLiteral("0"));
        xml.writeAttribute(QStringLiteral("file"), relativeFramePaths.at(index));
    }
    xml.writeEndElement();
    xml.writeEndElement();
    xml.writeEndDocument();
    if (!file.commit()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Failed to finalize PVD file: %1").arg(file.errorString());
        }
        return false;
    }
    return true;
}

streamcenterplus::io::ReadOptions canonicalNek5000ReadOptions() {
    streamcenterplus::io::ReadOptions options;
    // FlowData inspects the same canonical element-local GLL source that the
    // SEM evaluator and GUI use; no VTK cell topology is constructed while
    // metadata or fields are being probed.
    options.attachNek5000SpectralMetadata = true;
    return options;
}

bool enumerateNek5000Frames(const QString& filePath,
                            std::vector<streamcenterplus::io::DataSetFrame>* frames,
                            QString* errorMessage) {
    if (frames == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: Nek5000 frame storage was not provided.");
        }
        return false;
    }
    frames->clear();
    try {
        *frames = streamcenterplus::io::EnumerateDataSetFrames(fsPathFromQString(filePath));
    } catch (const std::exception& ex) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot enumerate Nek5000 time steps from %1:\n%2")
                                .arg(filePath, QString::fromLocal8Bit(ex.what()));
        }
        return false;
    }
    if (frames->empty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Nek5000 descriptor does not advertise any readable frames: %1")
                                .arg(filePath);
        }
        return false;
    }
    return true;
}

vtkSmartPointer<vtkDataSet> readCanonicalNek5000Frame(
    const streamcenterplus::io::DataSetFrame& frame,
    QString* errorMessage) {
    try {
        vtkSmartPointer<vtkDataSet> data = streamcenterplus::io::ReadDataSet(
            frame, canonicalNek5000ReadOptions());
        if (data == nullptr || !streamcenterplus::io::HasNek5000SpectralMetadata(data)) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("Nek5000 frame did not produce canonical spectral metadata: %1")
                                    .arg(qStringFromFsPath(frame.path));
            }
            return nullptr;
        }
        return data;
    } catch (const std::exception& ex) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot read canonical Nek5000 frame from %1:\n%2")
                                .arg(qStringFromFsPath(frame.path),
                                     QString::fromLocal8Bit(ex.what()));
        }
        return nullptr;
    }
}

bool canonicalNek5000Bounds(vtkDataSet* dataSet,
                            double bounds[6],
                            QString* errorMessage) {
    auto* pointSet = vtkPointSet::SafeDownCast(dataSet);
    vtkPoints* points = pointSet != nullptr ? pointSet->GetPoints() : nullptr;
    vtkDataArray* coordinates = points != nullptr ? points->GetData() : nullptr;
    if (coordinates == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "Canonical Nek5000 data does not expose a point-coordinate array.");
        }
        return false;
    }
    try {
        const streamcenterplus::io::Nek5000SpectralPointEvaluator geometry(
            dataSet, coordinates);
        // The evaluator publishes a conservative enclosure of the complete
        // high-order map.  It is intentionally safer than a raw GLL-node AABB
        // when a curved element reaches beyond its nodal coordinate extrema.
        const std::array<double, 6> spectralBounds = geometry.bounds();
        if (!std::all_of(spectralBounds.cbegin(), spectralBounds.cend(),
                         [](double value) { return std::isfinite(value); })) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                    "Canonical Nek5000 spectral geometry produced non-finite bounds.");
            }
            return false;
        }
        std::copy(spectralBounds.cbegin(), spectralBounds.cend(), bounds);
        return true;
    } catch (const std::exception& ex) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot evaluate canonical Nek5000 spectral bounds:\n%1")
                                .arg(QString::fromLocal8Bit(ex.what()));
        }
        return false;
    }
}

vtkSmartPointer<vtkDataSet> readAnyVtk(const QString& filePath, QString* errorMessage) {
    const QString ext = vtkExtensionSuffix(filePath);
    vtkDataSet* dataSet = nullptr;

    if (isTecplotExtension(ext)) {
        try {
            const std::vector<streamcenterplus::resample::FrameInfo> frames =
                streamcenterplus::resample::EnumerateDataSetFrames(
                    fsPathFromQString(filePath));
            if (frames.empty()) {
                throw std::runtime_error("The Tecplot source contains no frames.");
            }
            return streamcenterplus::resample::ReadDataSet(frames.front());
        } catch (const std::exception& ex) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("Cannot read Tecplot dataset %1:\n%2")
                                    .arg(filePath,
                                         QString::fromLocal8Bit(ex.what()));
            }
            return nullptr;
        }
    }

    try {
        if (ext == QStringLiteral("vtkhdf")) {
            vtkSmartPointer<vtkDataSet> data = streamcenterplus::vtkhdf::ReadDataSet(fsPathFromQString(filePath));
            if (data != nullptr) {
                return data;
            }
        }
    } catch (const std::exception& ex) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot read VTKHDF dataset %1:\n%2").arg(filePath, QString::fromLocal8Bit(ex.what()));
        }
        return nullptr;
    }

    if (isNek5000DescriptorExtension(ext)) {
        std::vector<streamcenterplus::io::DataSetFrame> frames;
        if (!enumerateNek5000Frames(filePath, &frames, errorMessage)) {
            return nullptr;
        }
        return readCanonicalNek5000Frame(frames.front(), errorMessage);
    }

    const QString realPath = realVtkPath(filePath);
    const QByteArray realPathBytes = QFile::encodeName(realPath);
    if (ext == QStringLiteral("vtk")) {
        auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(realPathBytes.constData());
        reader->Update();
        dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
        if (dataSet != nullptr) {
            vtkSmartPointer<vtkDataSet> retained;
            retained.TakeReference(vtkDataSet::SafeDownCast(dataSet->NewInstance()));
            retained->ShallowCopy(dataSet);
            return retained;
        }
    } else {
        auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
        reader->SetFileName(realPathBytes.constData());
        reader->Update();
        dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
        if (dataSet != nullptr) {
            vtkSmartPointer<vtkDataSet> retained;
            retained.TakeReference(vtkDataSet::SafeDownCast(dataSet->NewInstance()));
            retained->ShallowCopy(dataSet);
            return retained;
        }
    }

    if (errorMessage != nullptr) {
        *errorMessage = QString("Unsupported or unreadable VTK dataset: %1").arg(filePath);
    }
    return nullptr;
}

bool dimensionsForDataSet(vtkDataSet* dataSet, int dims[3]) {
    dims[0] = dims[1] = dims[2] = 0;
    if (auto* image = vtkImageData::SafeDownCast(dataSet)) {
        image->GetDimensions(dims);
        return true;
    }
    if (auto* rect = vtkRectilinearGrid::SafeDownCast(dataSet)) {
        rect->GetDimensions(dims);
        return true;
    }
    if (auto* structured = vtkStructuredGrid::SafeDownCast(dataSet)) {
        structured->GetDimensions(dims);
        return true;
    }
    return false;
}

bool axisIsRegular(const QVector<double>& values) {
    if (values.size() <= 1) {
        return true;
    }
    const double spacing = values[1] - values[0];
    if (!std::isfinite(spacing) || spacing == 0.0) {
        return false;
    }
    if (values.size() == 2) {
        return true;
    }
    for (int index = 2; index < values.size(); ++index) {
        const double expected = std::fma(static_cast<double>(index),
                                         spacing, values.first());
        const double magnitude = std::max(std::abs(values[index]),
                                          std::abs(expected));
        const double tolerance = std::abs(spacing)
                * static_cast<double>(index) * 1e-7
            + std::numeric_limits<double>::epsilon() * magnitude * 8.0
            + std::numeric_limits<double>::denorm_min() * 16.0;
        if (!std::isfinite(values[index]) || !std::isfinite(expected)
            || std::abs(values[index] - expected) > tolerance) {
            return false;
        }
    }
    return true;
}

QVector<double> vtkCoordinateValues(vtkDataArray* coordinates) {
    QVector<double> values;
    if (coordinates == nullptr || coordinates->GetNumberOfComponents() != 1
        || coordinates->GetNumberOfTuples() > std::numeric_limits<int>::max()) {
        return values;
    }
    values.reserve(static_cast<int>(coordinates->GetNumberOfTuples()));
    for (vtkIdType index = 0; index < coordinates->GetNumberOfTuples(); ++index) {
        const double value = coordinates->GetComponent(index, 0);
        if (!std::isfinite(value)) {
            return {};
        }
        values.append(value);
    }
    return values;
}

bool coordinatesEquivalent(double actual, double expected) {
    const double magnitude = std::max(std::abs(actual), std::abs(expected));
    const double tolerance = std::numeric_limits<double>::epsilon()
            * magnitude * 8.0
        + std::numeric_limits<double>::denorm_min() * 16.0;
    return std::isfinite(actual) && std::isfinite(expected)
        && std::abs(actual - expected) <= tolerance;
}

double ulpAtMagnitude(double magnitude) {
    magnitude = std::abs(magnitude);
    if (!std::isfinite(magnitude)) {
        return std::numeric_limits<double>::infinity();
    }
    if (magnitude == 0.0) {
        return std::numeric_limits<double>::denorm_min();
    }
    const double next = std::nextafter(
        magnitude, std::numeric_limits<double>::infinity());
    return std::isfinite(next)
        ? std::max(std::numeric_limits<double>::denorm_min(), next - magnitude)
        : std::numeric_limits<double>::epsilon() * magnitude;
}

double dataSetCoplanarityTolerance(vtkDataSet* dataSet) {
    double bounds[6] = {};
    dataSet->GetBounds(bounds);
    double geometricScale = 0.0;
    double coordinateMagnitude = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        geometricScale = std::max(
            geometricScale, std::abs(bounds[axis * 2 + 1] - bounds[axis * 2]));
        coordinateMagnitude = std::max(
            {coordinateMagnitude, std::abs(bounds[axis * 2]),
             std::abs(bounds[axis * 2 + 1])});
    }
    const double coordinateUlp = ulpAtMagnitude(coordinateMagnitude);
    return geometricScale * 1.0e-12 + coordinateUlp * 8.0;
}

bool dataSetPointsAreCoplanar(vtkDataSet* dataSet) {
    if (dataSet == nullptr || dataSet->GetNumberOfPoints() <= 3) {
        return true;
    }
    double bounds[6] = {};
    dataSet->GetBounds(bounds);
    for (int axis = 0; axis < 3; ++axis) {
        if (bounds[axis * 2] == bounds[axis * 2 + 1]) {
            return true;
        }
    }

    const double tolerance = dataSetCoplanarityTolerance(dataSet);
    double origin[3] = {};
    dataSet->GetPoint(0, origin);
    auto finitePoint = [](const double point[3]) {
        return std::isfinite(point[0]) && std::isfinite(point[1])
            && std::isfinite(point[2]);
    };
    if (!finitePoint(origin) || !std::isfinite(tolerance)) {
        return false;
    }

    vtkIdType farthestId = 0;
    double farthestDistance = 0.0;
    double point[3] = {};
    for (vtkIdType id = 1; id < dataSet->GetNumberOfPoints(); ++id) {
        dataSet->GetPoint(id, point);
        if (!finitePoint(point)) {
            return false;
        }
        const double distance = std::hypot(
            std::hypot(point[0] - origin[0], point[1] - origin[1]),
            point[2] - origin[2]);
        if (!std::isfinite(distance)) {
            return false;
        }
        if (distance > farthestDistance) {
            farthestDistance = distance;
            farthestId = id;
        }
    }
    if (farthestDistance <= tolerance) {
        return true;
    }

    double farthest[3] = {};
    dataSet->GetPoint(farthestId, farthest);
    const double direction[3] = {
        (farthest[0] - origin[0]) / farthestDistance,
        (farthest[1] - origin[1]) / farthestDistance,
        (farthest[2] - origin[2]) / farthestDistance,
    };
    double greatestLineDistance = 0.0;
    double planeRejection[3] = {};
    for (vtkIdType id = 1; id < dataSet->GetNumberOfPoints(); ++id) {
        dataSet->GetPoint(id, point);
        const double delta[3] = {point[0] - origin[0], point[1] - origin[1],
                                 point[2] - origin[2]};
        const double projection = delta[0] * direction[0]
            + delta[1] * direction[1] + delta[2] * direction[2];
        const double rejection[3] = {
            delta[0] - projection * direction[0],
            delta[1] - projection * direction[1],
            delta[2] - projection * direction[2],
        };
        const double distance = std::hypot(
            std::hypot(rejection[0], rejection[1]), rejection[2]);
        if (!std::isfinite(distance)) {
            return false;
        }
        if (distance > greatestLineDistance) {
            greatestLineDistance = distance;
            for (int axis = 0; axis < 3; ++axis) {
                planeRejection[axis] = rejection[axis];
            }
        }
    }
    if (greatestLineDistance <= tolerance) {
        return true;
    }

    const double secondDirection[3] = {
        planeRejection[0] / greatestLineDistance,
        planeRejection[1] / greatestLineDistance,
        planeRejection[2] / greatestLineDistance,
    };
    const double normal[3] = {
        direction[1] * secondDirection[2] - direction[2] * secondDirection[1],
        direction[2] * secondDirection[0] - direction[0] * secondDirection[2],
        direction[0] * secondDirection[1] - direction[1] * secondDirection[0],
    };
    for (vtkIdType id = 1; id < dataSet->GetNumberOfPoints(); ++id) {
        dataSet->GetPoint(id, point);
        const double distance = std::abs(
            (point[0] - origin[0]) * normal[0]
            + (point[1] - origin[1]) * normal[1]
            + (point[2] - origin[2]) * normal[2]);
        if (!std::isfinite(distance) || distance > tolerance) {
            return false;
        }
    }
    return true;
}

int spatialDimensionForDataSet(vtkDataSet* dataSet, const int dims[3]) {
    if (dataSet == nullptr) {
        return 0;
    }
    if (vtkImageData::SafeDownCast(dataSet) != nullptr
        || vtkRectilinearGrid::SafeDownCast(dataSet) != nullptr
        || vtkStructuredGrid::SafeDownCast(dataSet) != nullptr) {
        return dims[2] > 1 ? 3 : 2;
    }

    auto cellTypes = vtkSmartPointer<vtkCellTypes>::New();
    dataSet->GetCellTypes(cellTypes);
    for (vtkIdType index = 0; index < cellTypes->GetNumberOfTypes(); ++index) {
        if (vtkCellTypes::GetDimension(cellTypes->GetCellType(index)) >= 3) {
            return 3;
        }
    }
    return dataSetPointsAreCoplanar(dataSet) ? 2 : 3;
}

bool isUniformStructured(vtkDataSet* dataSet, int dims[3]) {
    if (auto* image = vtkImageData::SafeDownCast(dataSet)) {
        image->GetDimensions(dims);
        const double* spacing = image->GetSpacing();
        for (int axis = 0; axis < 3; ++axis) {
            if (dims[axis] > 1
                && (!std::isfinite(spacing[axis]) || spacing[axis] == 0.0)) {
                return false;
            }
        }
        return true;
    }
    if (auto* rectilinear = vtkRectilinearGrid::SafeDownCast(dataSet)) {
        rectilinear->GetDimensions(dims);
        const QVector<double> x = vtkCoordinateValues(rectilinear->GetXCoordinates());
        const QVector<double> y = vtkCoordinateValues(rectilinear->GetYCoordinates());
        const QVector<double> z = vtkCoordinateValues(rectilinear->GetZCoordinates());
        return x.size() == dims[0] && y.size() == dims[1] && z.size() == dims[2]
            && axisIsRegular(x) && axisIsRegular(y) && axisIsRegular(z);
    }
    auto* structured = vtkStructuredGrid::SafeDownCast(dataSet);
    if (structured == nullptr) {
        return false;
    }
    structured->GetDimensions(dims);
    if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) {
        return false;
    }
    QVector<double> axes[3];
    axes[0].reserve(dims[0]);
    axes[1].reserve(dims[1]);
    axes[2].reserve(dims[2]);
    double point[3] = {};
    for (int i = 0; i < dims[0]; ++i) {
        structured->GetPoint(i, point);
        axes[0].append(point[0]);
    }
    for (int j = 0; j < dims[1]; ++j) {
        structured->GetPoint(static_cast<vtkIdType>(j) * dims[0], point);
        axes[1].append(point[1]);
    }
    for (int k = 0; k < dims[2]; ++k) {
        structured->GetPoint(static_cast<vtkIdType>(k) * dims[0] * dims[1], point);
        axes[2].append(point[2]);
    }
    if (!axisIsRegular(axes[0]) || !axisIsRegular(axes[1])
        || !axisIsRegular(axes[2])) {
        return false;
    }
    for (int k = 0; k < dims[2]; ++k) {
        for (int j = 0; j < dims[1]; ++j) {
            for (int i = 0; i < dims[0]; ++i) {
                const vtkIdType id = i + static_cast<vtkIdType>(dims[0])
                    * (j + static_cast<vtkIdType>(dims[1]) * k);
                structured->GetPoint(id, point);
                if (!coordinatesEquivalent(point[0], axes[0].at(i))
                    || !coordinatesEquivalent(point[1], axes[1].at(j))
                    || !coordinatesEquivalent(point[2], axes[2].at(k))) {
                    return false;
                }
            }
        }
    }
    return true;
}

struct ArrayInfo {
    enum class Association {
        Unknown,
        Point,
        Cell,
    };

    QString name;
    int components = 0;
    vtkIdType tuples = 0;
    Association association = Association::Unknown;
};

void collectArrays(vtkFieldData* data,
                   ArrayInfo::Association association,
                   QList<ArrayInfo>* arrays,
                   QStringList* candidates) {
    if (data == nullptr || arrays == nullptr || candidates == nullptr) {
        return;
    }
    for (int index = 0; index < data->GetNumberOfArrays(); ++index) {
        vtkAbstractArray* abstractArray = data->GetAbstractArray(index);
        auto* array = vtkDataArray::SafeDownCast(abstractArray);
        if (array == nullptr || array->GetName() == nullptr) {
            continue;
        }
        ArrayInfo info;
        info.name = QString::fromUtf8(array->GetName()).trimmed();
        info.components = array->GetNumberOfComponents();
        info.tuples = array->GetNumberOfTuples();
        info.association = association;
        if (info.name.isEmpty()) {
            continue;
        }
        arrays->append(info);
        if (!candidates->contains(info.name, Qt::CaseInsensitive)) {
            candidates->append(info.name);
        }
        if (info.components >= 2) {
            const int componentCount = std::min(info.components, 3);
            for (int component = 0; component < componentCount; ++component) {
                const QString componentName = QStringLiteral("%1[%2]").arg(info.name).arg(component);
                if (!candidates->contains(componentName, Qt::CaseInsensitive)) {
                    candidates->append(componentName);
                }
            }
        }
    }
}

QStringList actualArrayNames(vtkFieldData* data) {
    QStringList names;
    if (data == nullptr) {
        return names;
    }
    for (int index = 0; index < data->GetNumberOfArrays(); ++index) {
        vtkAbstractArray* abstractArray = data->GetAbstractArray(index);
        auto* array = vtkDataArray::SafeDownCast(abstractArray);
        if (array == nullptr || array->GetName() == nullptr) {
            continue;
        }
        const QString name = QString::fromUtf8(array->GetName()).trimmed();
        if (!name.isEmpty() && !names.contains(name, Qt::CaseInsensitive)) {
            names.append(name);
        }
    }
    return names;
}

QStringList actualArrayNames(vtkDataSet* dataSet) {
    QStringList names;
    if (dataSet == nullptr) {
        return names;
    }
    for (const QString& name : actualArrayNames(dataSet->GetPointData())) {
        if (!names.contains(name, Qt::CaseInsensitive)) {
            names.append(name);
        }
    }
    for (const QString& name : actualArrayNames(dataSet->GetCellData())) {
        if (!names.contains(name, Qt::CaseInsensitive)) {
            names.append(name);
        }
    }
    names.sort(Qt::CaseInsensitive);
    return names;
}

QStringList intersectFieldNames(QStringList fields, const QStringList& next) {
    QStringList retained;
    for (const QString& field : std::as_const(fields)) {
        if (next.contains(field, Qt::CaseInsensitive) && !retained.contains(field, Qt::CaseInsensitive)) {
            retained.append(field);
        }
    }
    retained.sort(Qt::CaseInsensitive);
    return retained;
}

QList<ArrayInfo> filterArrayInfos(const QList<ArrayInfo>& arrays, const QStringList& allowedFields) {
    if (allowedFields.isEmpty()) {
        return arrays;
    }
    QList<ArrayInfo> filtered;
    for (const ArrayInfo& array : arrays) {
        if (allowedFields.contains(array.name, Qt::CaseInsensitive)) {
            filtered.append(array);
        }
    }
    return filtered;
}

void addCandidatesForArray(const ArrayInfo& info, QStringList* candidates) {
    if (candidates == nullptr || info.name.trimmed().isEmpty()) {
        return;
    }
    if (!candidates->contains(info.name, Qt::CaseInsensitive)) {
        candidates->append(info.name);
    }
    if (info.components >= 2) {
        const int componentCount = std::min(info.components, 3);
        for (int component = 0; component < componentCount; ++component) {
            const QString componentName = QStringLiteral("%1[%2]").arg(info.name).arg(component);
            if (!candidates->contains(componentName, Qt::CaseInsensitive)) {
                candidates->append(componentName);
            }
        }
    }
}

QStringList candidateListForArrays(const QList<ArrayInfo>& arrays) {
    QStringList candidates;
    for (const ArrayInfo& array : arrays) {
        addCandidatesForArray(array, &candidates);
    }
    candidates.sort(Qt::CaseInsensitive);
    return candidates;
}

bool sameCandidateFieldSet(const QStringList& lhs, const QStringList& rhs) {
    return lhs.size() == rhs.size()
        && std::all_of(lhs.cbegin(), lhs.cend(), [&](const QString& candidate) {
               return rhs.contains(candidate, Qt::CaseSensitive);
           });
}

bool sameArraySchema(const QList<ArrayInfo>& lhs,
                     const QList<ArrayInfo>& rhs,
                     bool requireSameTupleCount) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    QVector<bool> matched(rhs.size(), false);
    for (const ArrayInfo& left : lhs) {
        bool found = false;
        for (qsizetype index = 0; index < rhs.size(); ++index) {
            const ArrayInfo& right = rhs.at(index);
            if (!matched.at(index)
                && left.name.compare(right.name, Qt::CaseSensitive) == 0
                && left.components == right.components
                && (!requireSameTupleCount || left.tuples == right.tuples)
                && left.association == right.association) {
                matched[index] = true;
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

struct VtkDataSummary {
    bool hasBounds = false;
    double bounds[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int dims[3] = {0, 0, 0};
    bool hasDims = false;
    bool uniformStructured = false;
    int spatialDimension = 0;
    QList<ArrayInfo> arrays;
    QStringList termCandidates;
    QList<QByteArray> coordinateSignatures;
};

template <typename Value>
void addGeometryHashValue(QCryptographicHash* hash, const Value& value) {
    hash->addData(reinterpret_cast<const char*>(&value),
                  static_cast<qsizetype>(sizeof(value)));
}

QByteArray coordinateSignature(vtkDataSet* dataSet) {
    QCryptographicHash hash(QCryptographicHash::Sha256);
    constexpr char version[] = "SCP-FlowDataCoordinates-v1";
    hash.addData(version, static_cast<qsizetype>(sizeof(version) - 1));
    if (dataSet == nullptr) {
        return hash.result();
    }

    const vtkIdType pointCount = dataSet->GetNumberOfPoints();
    std::vector<std::array<double, 3>> points;
    points.reserve(static_cast<std::size_t>(pointCount));
    bool hasNonFiniteCoordinate = false;
    for (vtkIdType pointId = 0; pointId < pointCount; ++pointId) {
        std::array<double, 3> point{};
        dataSet->GetPoint(pointId, point.data());
        for (double& coordinate : point) {
            // Treat the two IEEE representations of physical zero equally.
            if (coordinate == 0.0) {
                coordinate = 0.0;
            }
            hasNonFiniteCoordinate = hasNonFiniteCoordinate || !std::isfinite(coordinate);
        }
        points.push_back(point);
    }
    if (hasNonFiniteCoordinate) {
        // Invalid geometry must never be made equal by an unordered NaN sort.
        // Preserve the source tuple sequence in this diagnostic signature.
        addGeometryHashValue(&hash, static_cast<qsizetype>(points.size()));
        for (const std::array<double, 3>& point : points) {
            hash.addData(reinterpret_cast<const char*>(point.data()),
                         static_cast<qsizetype>(sizeof(point)));
        }
        return hash.result();
    }
    std::sort(points.begin(), points.end(), [](const auto& lhs, const auto& rhs) {
        return std::lexicographical_compare(lhs.cbegin(), lhs.cend(),
                                            rhs.cbegin(), rhs.cend());
    });
    points.erase(std::unique(points.begin(), points.end()), points.end());
    addGeometryHashValue(&hash, static_cast<qsizetype>(points.size()));
    for (const std::array<double, 3>& point : points) {
        hash.addData(reinterpret_cast<const char*>(point.data()),
                     static_cast<qsizetype>(sizeof(point)));
    }
    return hash.result();
}

void mergeBounds(VtkDataSummary* summary, const double* bounds) {
    if (summary == nullptr || bounds == nullptr) {
        return;
    }
    if (!summary->hasBounds) {
        for (int index = 0; index < 6; ++index) {
            summary->bounds[index] = bounds[index];
        }
        summary->hasBounds = true;
        return;
    }
    summary->bounds[0] = std::min(summary->bounds[0], bounds[0]);
    summary->bounds[1] = std::max(summary->bounds[1], bounds[1]);
    summary->bounds[2] = std::min(summary->bounds[2], bounds[2]);
    summary->bounds[3] = std::max(summary->bounds[3], bounds[3]);
    summary->bounds[4] = std::min(summary->bounds[4], bounds[4]);
    summary->bounds[5] = std::max(summary->bounds[5], bounds[5]);
}

bool summarizeDataSet(vtkDataSet* dataSet, VtkDataSummary* summary) {
    const bool canonicalNek5000 = dataSet != nullptr
        && streamcenterplus::io::HasNek5000SpectralMetadata(dataSet);
    if (dataSet == nullptr || summary == nullptr
        || dataSet->GetNumberOfPoints() <= 0
        // Canonical Native Nek5000 represents each spectral element as a
        // local GLL point block.  It intentionally has no linear VTK cells;
        // rejecting it here would force FlowData probing onto an invalid
        // linearized representation.
        || (dataSet->GetNumberOfCells() <= 0 && !canonicalNek5000)) {
        return false;
    }

    double bounds[6] = {};
    dataSet->GetBounds(bounds);
    mergeBounds(summary, bounds);

    int dims[3] = {0, 0, 0};
    const bool hasDims = dimensionsForDataSet(dataSet, dims);
    if (!summary->hasDims && hasDims) {
        for (int index = 0; index < 3; ++index) {
            summary->dims[index] = dims[index];
        }
        summary->hasDims = true;
    }
    summary->spatialDimension = std::max(
        summary->spatialDimension, spatialDimensionForDataSet(dataSet, dims));
    summary->uniformStructured = summary->uniformStructured || isUniformStructured(dataSet, dims);
    summary->coordinateSignatures.append(coordinateSignature(dataSet));

    if (summary->arrays.isEmpty()) {
        collectArrays(dataSet->GetPointData(), ArrayInfo::Association::Point,
                      &summary->arrays, &summary->termCandidates);
        collectArrays(dataSet->GetCellData(), ArrayInfo::Association::Cell,
                      &summary->arrays, &summary->termCandidates);
        summary->termCandidates.sort(Qt::CaseInsensitive);
    }
    return true;
}

QStringList collectVtmBlockFiles(const QString& vtmPath, QString* errorMessage) {
    QStringList files;
    QFile file(vtmPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot read VTM file: %1").arg(vtmPath);
        }
        return files;
    }

    QXmlStreamReader xml(&file);
    const QDir baseDir(QFileInfo(vtmPath).absolutePath());
    while (!xml.atEnd()) {
        xml.readNext();
        if (!xml.isStartElement() || xml.name() != QStringLiteral("DataSet")) {
            continue;
        }
        const QString relativeFile = xml.attributes().value(QStringLiteral("file")).toString().trimmed();
        if (relativeFile.isEmpty()) {
            continue;
        }
        const QString path = QDir::isAbsolutePath(relativeFile)
            ? QFileInfo(relativeFile).absoluteFilePath()
            : baseDir.absoluteFilePath(relativeFile);
        files.append(normalizedPath(path));
    }
    if (xml.hasError()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot parse VTM file %1: %2").arg(vtmPath, xml.errorString());
        }
        return {};
    }
    if (files.isEmpty() && errorMessage != nullptr) {
        *errorMessage = QString("VTM file does not reference any VTK block: %1").arg(vtmPath);
    }
    return files;
}

bool coordinateBoundsEquivalent(const double* lhs, const double* rhs) {
    if (lhs == nullptr || rhs == nullptr) {
        return false;
    }
    for (int axis = 0; axis < 3; ++axis) {
        const int lower = axis * 2;
        const int upper = lower + 1;
        const double extent = std::max(
            std::abs(lhs[upper] - lhs[lower]),
            std::abs(rhs[upper] - rhs[lower]));
        const double coordinateMagnitude = std::max(
            {std::abs(lhs[lower]), std::abs(lhs[upper]),
             std::abs(rhs[lower]), std::abs(rhs[upper])});
        const double tolerance = extent * 1.0e-8
            + ulpAtMagnitude(coordinateMagnitude) * 8.0;
        for (const int index : {lower, upper}) {
            if (!std::isfinite(lhs[index]) || !std::isfinite(rhs[index])
                || std::abs(lhs[index] - rhs[index]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

bool coordinatesEquivalent(const VtkDataSummary& lhs, const VtkDataSummary& rhs) {
    if (!lhs.hasBounds || !rhs.hasBounds) {
        return false;
    }
    if (!coordinateBoundsEquivalent(lhs.bounds, rhs.bounds)) {
        return false;
    }
    if (lhs.hasDims != rhs.hasDims) {
        return false;
    }
    if (lhs.hasDims) {
        for (int index = 0; index < 3; ++index) {
            if (lhs.dims[index] != rhs.dims[index]) {
                return false;
            }
        }
    }
    return lhs.coordinateSignatures == rhs.coordinateSignatures;
}

bool summarizeAnyVtk(const QString& filePath, VtkDataSummary* summary, QString* errorMessage) {
    if (summary == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: VTK summary storage was not provided.");
        }
        return false;
    }
    *summary = VtkDataSummary{};

    const QString ext = vtkExtensionSuffix(filePath);
    if (ext == QStringLiteral("vtm")) {
        QString blockError;
        const QStringList blockFiles = collectVtmBlockFiles(filePath, &blockError);
        if (blockFiles.isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = blockError;
            }
            return false;
        }
        for (const QString& blockFile : blockFiles) {
            QString readError;
            vtkSmartPointer<vtkDataSet> dataSet = readAnyVtk(blockFile, &readError);
            if (dataSet == nullptr) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("Cannot read VTM block %1 referenced by %2:\n%3")
                                        .arg(blockFile, filePath, readError);
                }
                return false;
            }
            summarizeDataSet(dataSet, summary);
        }
        if (!summary->hasBounds || summary->termCandidates.isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("VTM file does not contain readable volume blocks with numeric arrays: %1").arg(filePath);
            }
            return false;
        }
        summary->uniformStructured = false;
        return true;
    }

    QString readError;
    vtkSmartPointer<vtkDataSet> dataSet = readAnyVtk(filePath, &readError);
    if (dataSet == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = readError;
        }
        return false;
    }
    if (!summarizeDataSet(dataSet, summary)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Unsupported or empty VTK dataset: %1").arg(filePath);
        }
        return false;
    }
    return true;
}

QString baseFieldName(QString field) {
    field = field.trimmed();
    const QRegularExpression componentRegex(QStringLiteral(R"(^(.*)\[\d+\]$)"));
    const QRegularExpressionMatch match = componentRegex.match(field);
    if (match.hasMatch()) {
        return match.captured(1).trimmed();
    }
    return field;
}

QStringList actualFieldsFromCandidates(const QStringList& fields) {
    QStringList names;
    for (const QString& field : fields) {
        const QString name = baseFieldName(field);
        if (!name.isEmpty() && !names.contains(name, Qt::CaseInsensitive)) {
            names.append(name);
        }
    }
    names.sort(Qt::CaseInsensitive);
    return names;
}

void filterFieldDataArrays(vtkFieldData* data, const QStringList& allowedFields) {
    if (data == nullptr || allowedFields.isEmpty()) {
        return;
    }
    for (int index = data->GetNumberOfArrays() - 1; index >= 0; --index) {
        vtkAbstractArray* array = data->GetAbstractArray(index);
        const char* rawName = array != nullptr ? array->GetName() : nullptr;
        const QString name = QString::fromUtf8(rawName == nullptr ? "" : rawName).trimmed();
        if (name.isEmpty() || !allowedFields.contains(name, Qt::CaseInsensitive)) {
            data->RemoveArray(index);
        }
    }
}

void filterDataSetArrays(vtkDataSet* dataSet, const QStringList& allowedFields) {
    if (dataSet == nullptr || allowedFields.isEmpty()) {
        return;
    }
    filterFieldDataArrays(dataSet->GetPointData(), allowedFields);
    filterFieldDataArrays(dataSet->GetCellData(), allowedFields);
}

const ArrayInfo* findArray(const QList<ArrayInfo>& arrays, const QString& name) {
    for (const ArrayInfo& array : arrays) {
        if (array.name.compare(name, Qt::CaseInsensitive) == 0) {
            return &array;
        }
    }
    return nullptr;
}

void inferTerms(const QList<ArrayInfo>& arrays, int dimension, FlowDataEntry* entry) {
    if (entry == nullptr) {
        return;
    }

    const QList<std::tuple<QString, QString, QString>> triples = {
        {QStringLiteral("u"), QStringLiteral("v"), QStringLiteral("w")},
        {QStringLiteral("U"), QStringLiteral("V"), QStringLiteral("W")},
        {QStringLiteral("Ux"), QStringLiteral("Uy"), QStringLiteral("Uz")},
        {QStringLiteral("Vx"), QStringLiteral("Vy"), QStringLiteral("Vz")},
        {QStringLiteral("mean_u"), QStringLiteral("mean_v"), QStringLiteral("mean_w")},
        {QStringLiteral("u_mean"), QStringLiteral("v_mean"), QStringLiteral("w_mean")},
    };
    for (const auto& triple : triples) {
        if (findArray(arrays, std::get<0>(triple)) != nullptr
            && findArray(arrays, std::get<1>(triple)) != nullptr
            && findArray(arrays, std::get<2>(triple)) != nullptr) {
            entry->uTerm = std::get<0>(triple);
            entry->vTerm = std::get<1>(triple);
            entry->wTerm = std::get<2>(triple);
            break;
        }
    }

    if (entry->uTerm.isEmpty() || entry->vTerm.isEmpty()) {
        const QList<QPair<QString, QString>> pairs = {
            {QStringLiteral("u"), QStringLiteral("v")},
            {QStringLiteral("U"), QStringLiteral("V")},
            {QStringLiteral("Ux"), QStringLiteral("Uy")},
            {QStringLiteral("Vx"), QStringLiteral("Vy")},
            {QStringLiteral("mean_u"), QStringLiteral("mean_v")},
            {QStringLiteral("u_mean"), QStringLiteral("v_mean")},
        };
        for (const auto& pair : pairs) {
            if (findArray(arrays, pair.first) != nullptr && findArray(arrays, pair.second) != nullptr) {
                entry->uTerm = pair.first;
                entry->vTerm = pair.second;
                entry->wTerm = dimension <= 2 ? QStringLiteral("0") : QString();
                break;
            }
        }
    }

    if (entry->uTerm.isEmpty() || entry->vTerm.isEmpty()) {
        const QStringList vectorCandidates = {
            QStringLiteral("velocity"),
            QStringLiteral("mean_velocity"),
            QStringLiteral("phi_velocity"),
            QStringLiteral("phi_velocity_real"),
            QStringLiteral("velocity_rom"),
            QStringLiteral("U"),
            QStringLiteral("V"),
            QStringLiteral("Velocity"),
            QStringLiteral("Vel"),
            QStringLiteral("vel"),
            QStringLiteral("data"),
            QStringLiteral("Data"),
        };
        for (const QString& name : vectorCandidates) {
            const ArrayInfo* array = findArray(arrays, name);
            if (array != nullptr && array->components >= 2) {
                entry->uTerm = QStringLiteral("%1[0]").arg(array->name);
                entry->vTerm = QStringLiteral("%1[1]").arg(array->name);
                entry->wTerm = array->components >= 3 ? QStringLiteral("%1[2]").arg(array->name) : QStringLiteral("0");
                break;
            }
        }
    }

    const QStringList pressureCandidates = {
        QStringLiteral("pressure"),
        QStringLiteral("p"),
        QStringLiteral("mean_pressure"),
        QStringLiteral("p_mean"),
        QStringLiteral("pi"),
        QStringLiteral("kinematic_pressure"),
    };
    for (const QString& name : pressureCandidates) {
        if (findArray(arrays, name) != nullptr) {
            entry->pTerm = name;
            break;
        }
    }
}

bool coordinatesEquivalent(vtkDataSet* lhs, vtkDataSet* rhs) {
    if (lhs == nullptr || rhs == nullptr) {
        return false;
    }
    double left[6] = {};
    double right[6] = {};
    lhs->GetBounds(left);
    rhs->GetBounds(right);
    if (!coordinateBoundsEquivalent(left, right)) {
        return false;
    }
    int leftDims[3] = {};
    int rightDims[3] = {};
    const bool leftHasDims = dimensionsForDataSet(lhs, leftDims);
    const bool rightHasDims = dimensionsForDataSet(rhs, rightDims);
    if (leftHasDims != rightHasDims) {
        return false;
    }
    if (leftHasDims) {
        for (int index = 0; index < 3; ++index) {
            if (leftDims[index] != rightDims[index]) {
                return false;
            }
        }
    }
    return coordinateSignature(lhs) == coordinateSignature(rhs);
}

QString boundsSummary(const double* bounds) {
    return QStringLiteral("X%1 Y%2 Z%3")
        .arg(formatRange(bounds[0], bounds[1]),
             formatRange(bounds[2], bounds[3]),
             formatRange(bounds[4], bounds[5]));
}

QString boundsSummary(vtkDataSet* dataSet) {
    double bounds[6] = {};
    dataSet->GetBounds(bounds);
    return boundsSummary(bounds);
}

Qt::CaseSensitivity fileSystemPathCaseSensitivity() {
#ifdef Q_OS_WIN
    return Qt::CaseInsensitive;
#else
    return Qt::CaseSensitive;
#endif
}

QString sourcePathIdentity(const QString& path) {
    const QFileInfo info(path);
    const QString canonical = info.canonicalFilePath();
    return normalizedPath(QDir::cleanPath(
        canonical.isEmpty() ? info.absoluteFilePath() : canonical));
}

bool sameSourcePath(const QString& lhs, const QString& rhs) {
    return sourcePathIdentity(lhs).compare(sourcePathIdentity(rhs),
                                            fileSystemPathCaseSensitivity()) == 0;
}

bool containsSourcePath(const QStringList& paths, const QString& candidate) {
    return std::any_of(paths.cbegin(), paths.cend(), [&](const QString& path) {
        return sameSourcePath(path, candidate);
    });
}

QString sourcePathIdentityKey(const QString& path) {
    const QString identity = sourcePathIdentity(path);
    return fileSystemPathCaseSensitivity() == Qt::CaseInsensitive
        ? identity.toCaseFolded()
        : identity;
}

bool appendBoundedCasePath(const QString& path,
                           QStringList* caseFiles,
                           QSet<QString>* identities,
                           qint64* aggregateBytes,
                           QString* errorMessage) {
    const QString identity = sourcePathIdentityKey(path);
    if (identities->contains(identity)) {
        return true;
    }
    const qsizetype pathBytes = path.toUtf8().size();
    if (pathBytes > Streamcenter::FlowDataProtocolLimits::MaximumStringBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("An EnSight Gold case path exceeds the %1-byte path limit: %2")
                                .arg(Streamcenter::FlowDataProtocolLimits::MaximumStringBytes)
                                .arg(path);
        }
        return false;
    }
    if (caseFiles->size() >= Streamcenter::FlowDataProtocolLimits::MaximumCaseCount) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold selection exceeds the %1-case limit.")
                                .arg(Streamcenter::FlowDataProtocolLimits::MaximumCaseCount);
        }
        return false;
    }
    if (*aggregateBytes + pathBytes
        > Streamcenter::FlowDataProtocolLimits::MaximumAggregateStringBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case paths exceed the %1-byte aggregate limit.")
                                .arg(Streamcenter::FlowDataProtocolLimits::MaximumAggregateStringBytes);
        }
        return false;
    }
    caseFiles->append(path);
    identities->insert(identity);
    *aggregateBytes += pathBytes;
    return true;
}

QStringList collectEnSightGoldCaseFiles(const QStringList& sources, QString* errorMessage) {
    QStringList caseFiles;
    QSet<QString> caseFileIdentities;
    qint64 aggregatePathBytes = 0;
    if (sources.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Select at least one EnSight Gold *.case file or a folder containing *.case files.");
        }
        return caseFiles;
    }

    for (const QString& source : sources) {
        if (source.isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("An EnSight Gold source path is empty.");
            }
            caseFiles.clear();
            return caseFiles;
        }
        const QFileInfo info(source);
        if (!info.exists()) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold source does not exist: %1").arg(source);
            }
            caseFiles.clear();
            return caseFiles;
        }

        if (info.isDir()) {
            QStringList directoryCaseFiles;
            QSet<QString> directoryIdentities;
            bool foundCaseEntry = false;
            QDirIterator iterator(info.absoluteFilePath(),
                                  {QStringLiteral("*.case")},
                                  QDir::Files | QDir::NoDotAndDotDot | QDir::Readable,
                                  QDirIterator::NoIteratorFlags);
            while (iterator.hasNext()) {
                const QString path = normalizedPath(QFileInfo(iterator.next()).absoluteFilePath());
                foundCaseEntry = true;
                const QString identity = sourcePathIdentityKey(path);
                if (caseFileIdentities.contains(identity)
                    || directoryIdentities.contains(identity)) {
                    continue;
                }
                if (caseFiles.size() + directoryCaseFiles.size()
                        >= Streamcenter::FlowDataProtocolLimits::MaximumCaseCount) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight Gold selection exceeds the %1-case limit.")
                                            .arg(Streamcenter::FlowDataProtocolLimits::MaximumCaseCount);
                    }
                    caseFiles.clear();
                    return caseFiles;
                }
                const qsizetype pathBytes = path.toUtf8().size();
                if (pathBytes > Streamcenter::FlowDataProtocolLimits::MaximumStringBytes
                    || aggregatePathBytes + pathBytes
                        > Streamcenter::FlowDataProtocolLimits::MaximumAggregateStringBytes) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight Gold case paths exceed the bounded protocol size.");
                    }
                    caseFiles.clear();
                    return caseFiles;
                }
                directoryCaseFiles.append(path);
                directoryIdentities.insert(identity);
                aggregatePathBytes += pathBytes;
            }
            if (!foundCaseEntry) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("No EnSight Gold case files (*.case) were found in folder: %1")
                                        .arg(info.absoluteFilePath());
                }
                caseFiles.clear();
                return caseFiles;
            }
            directoryCaseFiles.sort(fileSystemPathCaseSensitivity());
            caseFiles.append(directoryCaseFiles);
            caseFileIdentities.unite(directoryIdentities);
            continue;
        }

        if (!info.isFile() || info.suffix().compare(QStringLiteral("case"), Qt::CaseInsensitive) != 0) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold import expects *.case files or folders containing *.case files: %1")
                                    .arg(source);
            }
            caseFiles.clear();
            return caseFiles;
        }

        const QString path = normalizedPath(info.absoluteFilePath());
        if (!appendBoundedCasePath(path,
                                   &caseFiles,
                                   &caseFileIdentities,
                                   &aggregatePathBytes,
                                   errorMessage)) {
            caseFiles.clear();
            return caseFiles;
        }
    }

    return caseFiles;
}

struct EnSightGoldTimeSet {
    int id = 0;
    int declarationLine = 0;
    int declaredSteps = 0;
    bool hasDeclaredSteps = false;
    bool hasExplicitFilenameNumbers = false;
    bool hasFilenameStart = false;
    bool hasFilenameIncrement = false;
    qint64 filenameStart = 0;
    qint64 filenameIncrement = 0;
    QVector<qint64> filenameNumbers;
    QVector<double> timeValues;
};

struct EnSightGoldFileChunk {
    bool hasFilenameIndex = false;
    qint64 filenameIndex = 0;
    int stepCount = 0;
};

struct EnSightGoldFileSet {
    int id = 0;
    int declarationLine = 0;
    bool usesFilenameIndices = false;
    QVector<EnSightGoldFileChunk> chunks;
};

struct EnSightGoldReference {
    enum class Role {
        ModelGeometry,
        MeasuredGeometry,
        VolumeVariable,
        MeasuredVariable,
        ConstantVariable,
        Auxiliary,
    };

    Role role = Role::Auxiliary;
    QString declaration;
    int declarationLine = 0;
    int timeSetId = 0;
    int fileSetId = 0;
    QString variableName;
    bool hasFrequency = false;
    bool frequencyUndefined = false;
    double frequency = 0.0;
    QStringList patterns;
    QStringList resolvedPaths;
    QVector<QStringList> resolvedPatternPaths;
    QVector<QStringList> expandedPatternNames;
};

struct EnSightGoldConstantPerCase {
    int declarationLine = 0;
    int timeSetId = 0;
    QString variableName;
    QVector<double> values;
};

struct EnSightGoldCaseSource {
    QString path;
    QVector<double> timeValues;
    QList<ArrayInfo> arrays;
    QHash<int, EnSightGoldTimeSet> timeSets;
    QHash<int, EnSightGoldFileSet> fileSets;
    QVector<EnSightGoldReference> references;
    QVector<EnSightGoldConstantPerCase> constantsPerCase;
    QStringList constantPerCaseDeclarations;
    QStringList referencedPaths;
    QStringList modelGeometryPaths;
    int maximumTimeSteps = 0;
    bool hasMaximumTimeSteps = false;
    int blockContinuationCount = 0;
    qint64 blockContinuationStart = 0;
    qint64 blockContinuationIncrement = 0;
    bool hasBlockContinuationCount = false;
    bool hasBlockContinuationStart = false;
    bool hasBlockContinuationIncrement = false;
};

struct EnSightGoldFrameSource {
    int caseIndex = -1;
    double timeValue = 0.0;
};

constexpr qint64 kLargeEnSightGoldSingleFileBytes = 512ll * 1024ll * 1024ll;
constexpr qint64 kLargeEnSightGoldFrameBytes = 1024ll * 1024ll * 1024ll;

QString cleanEnSightGoldName(QString name) {
    name = name.trimmed();
    if (name.size() >= 2
        && ((name.front() == QLatin1Char('"') && name.back() == QLatin1Char('"'))
            || (name.front() == QLatin1Char('\'') && name.back() == QLatin1Char('\'')))) {
        name = name.mid(1, name.size() - 2).trimmed();
    }
    return name;
}

bool normalizeEnSightGoldVariableName(const QString& rawName,
                                      QString* normalizedName,
                                      QString* errorMessage,
                                      int lineNumber) {
    if (normalizedName == nullptr) {
        return false;
    }
    QString name = cleanEnSightGoldName(rawName);
    name.replace(QRegularExpression(QStringLiteral(R"(\s+)")), QStringLiteral("_"));
    static const QRegularExpression safeName(
        QStringLiteral(R"(^[\p{L}_][\p{L}\p{N}_]*$)"),
        QRegularExpression::UseUnicodePropertiesOption);
    if (name.size() > 79 || !safeName.match(name).hasMatch()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold variable name on line %1 must be at most 79 characters, must not start with a digit, and must not contain whitespace/reserved characters; quoted spaces are normalized to underscores.")
                                .arg(lineNumber);
        }
        return false;
    }
    *normalizedName = name;
    return true;
}

bool parseIntToken(QString text, int* value) {
    if (value == nullptr) {
        return false;
    }
    text = text.trimmed();
    bool ok = false;
    const int parsed = text.toInt(&ok);
    if (ok) {
        *value = parsed;
    }
    return ok;
}

bool parsePositiveIntToken(const QString& text, int* value) {
    int parsed = 0;
    if (!parseIntToken(text, &parsed) || parsed <= 0) {
        return false;
    }
    if (value != nullptr) {
        *value = parsed;
    }
    return true;
}

bool parseEnSightGoldIntegerToken(const QString& text, qint64* value) {
    if (value == nullptr) {
        return false;
    }
    bool ok = false;
    const qint64 parsed = text.trimmed().toLongLong(&ok);
    if (!ok || parsed < std::numeric_limits<int>::min()
        || parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    *value = parsed;
    return true;
}

struct EnSightGoldCaseLine {
    int number = 0;
    QString text;
};

QString stripEnSightGoldInlineComment(const QString& rawLine) {
    QChar quote;
    for (qsizetype index = 0; index < rawLine.size(); ++index) {
        const QChar character = rawLine.at(index);
        if (!quote.isNull()) {
            if (character == quote) {
                quote = QChar{};
            }
            continue;
        }
        if (character == QLatin1Char('"') || character == QLatin1Char('\'')) {
            quote = character;
        } else if (character == QLatin1Char('#')) {
            return rawLine.left(index);
        }
    }
    return rawLine;
}

bool readBoundedEnSightGoldCaseLines(const QByteArray& payload,
                                     QVector<EnSightGoldCaseLine>* lines,
                                     QString* errorMessage) {
    if (lines == nullptr) {
        return false;
    }
    lines->clear();
    constexpr qsizetype kMaximumCaseLineBytes = 1024;
    qsizetype offset = 0;
    int lineNumber = 1;
    while (offset <= payload.size()) {
        qsizetype end = payload.indexOf('\n', offset);
        if (end < 0) {
            end = payload.size();
        }
        qsizetype length = end - offset;
        if (length > 0 && payload.at(offset + length - 1) == '\r') {
            --length;
        }
        if (length > kMaximumCaseLineBytes) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold case line %1 exceeds the 1024-byte format limit.")
                                    .arg(lineNumber);
            }
            return false;
        }
        const QByteArray rawBytes = payload.mid(offset, length);
        if (rawBytes.contains('\0')) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold case line %1 contains a NUL byte.")
                                    .arg(lineNumber);
            }
            return false;
        }
        QString text = stripEnSightGoldInlineComment(QString::fromUtf8(rawBytes)).trimmed();
        if (!text.isEmpty()) {
            lines->append(EnSightGoldCaseLine{lineNumber, text});
        }
        if (end == payload.size()) {
            break;
        }
        offset = end + 1;
        ++lineNumber;
    }
    return true;
}

bool tokenizeEnSightGoldPayload(const QString& payload,
                                QStringList* tokens,
                                QString* errorMessage,
                                int lineNumber) {
    if (tokens == nullptr) {
        return false;
    }
    tokens->clear();
    QString token;
    QChar quote;
    bool tokenStarted = false;
    for (qsizetype index = 0; index < payload.size(); ++index) {
        const QChar character = payload.at(index);
        if (!quote.isNull()) {
            if (character == quote) {
                quote = QChar{};
            } else {
                token.append(character);
            }
            tokenStarted = true;
            continue;
        }
        if (character == QLatin1Char('"') || character == QLatin1Char('\'')) {
            if (tokenStarted) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "EnSight Gold case line %1 contains a quote in the middle of a token.")
                                        .arg(lineNumber);
                }
                return false;
            }
            quote = character;
            tokenStarted = true;
        } else if (character.isSpace() || character == QLatin1Char(',')) {
            if (tokenStarted) {
                tokens->append(token);
                token.clear();
                tokenStarted = false;
            }
        } else {
            token.append(character);
            tokenStarted = true;
        }
    }
    if (!quote.isNull()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case line %1 contains an unterminated quoted token.")
                                .arg(lineNumber);
        }
        return false;
    }
    if (tokenStarted) {
        tokens->append(token);
    }
    return true;
}

int enSightGoldUnquotedColon(const QString& line) {
    QChar quote;
    for (qsizetype index = 0; index < line.size(); ++index) {
        const QChar character = line.at(index);
        if (!quote.isNull()) {
            if (character == quote) {
                quote = QChar{};
            }
            continue;
        }
        if (character == QLatin1Char('"') || character == QLatin1Char('\'')) {
            quote = character;
        } else if (character == QLatin1Char(':')) {
            return static_cast<int>(index);
        }
    }
    return -1;
}

int enSightGoldComponentCount(const QString& declaration) {
    const QString lower = declaration.toLower();
    if (lower.contains(QStringLiteral("tensor asym"))) {
        return 9;
    }
    if (lower.contains(QStringLiteral("tensor symm")) || lower.contains(QStringLiteral("tensor"))) {
        return 6;
    }
    if (lower.contains(QStringLiteral("vector"))) {
        return 3;
    }
    if (lower.contains(QStringLiteral("complex"))) {
        return 2;
    }
    return 1;
}

QString enSightGoldSectionName(const QString& line) {
    if (line == QStringLiteral("FORMAT")
        || line == QStringLiteral("GEOMETRY")
        || line == QStringLiteral("VARIABLE")
        || line == QStringLiteral("TIME")
        || line == QStringLiteral("FILE")
        || line == QStringLiteral("SCRIPTS")
        || line == QStringLiteral("BLOCK_CONTINUATION")
        || line == QStringLiteral("QUERY")
        || line == QStringLiteral("MATERIAL")) {
        return line;
    }
    return {};
}

bool isEnSightGoldSectionSpelling(const QString& line) {
    const QString upper = line.toUpper();
    return upper == QStringLiteral("FORMAT")
        || upper == QStringLiteral("GEOMETRY")
        || upper == QStringLiteral("VARIABLE")
        || upper == QStringLiteral("TIME")
        || upper == QStringLiteral("FILE")
        || upper == QStringLiteral("SCRIPTS")
        || upper == QStringLiteral("BLOCK_CONTINUATION")
        || upper == QStringLiteral("QUERY")
        || upper == QStringLiteral("MATERIAL");
}

bool takeEnSightGoldBinding(QStringList* tokens,
                            qsizetype minimumTailTokens,
                            bool allowFileSet,
                            int* timeSetId,
                            int* fileSetId) {
    if (tokens == nullptr || timeSetId == nullptr || fileSetId == nullptr) {
        return false;
    }
    *timeSetId = 0;
    *fileSetId = 0;
    int parsed = 0;
    if (tokens->size() > minimumTailTokens
        && parsePositiveIntToken(tokens->first(), &parsed)) {
        *timeSetId = parsed;
        tokens->removeFirst();
        if (allowFileSet && tokens->size() > minimumTailTokens
            && parsePositiveIntToken(tokens->first(), &parsed)) {
            *fileSetId = parsed;
            tokens->removeFirst();
        }
    }
    return true;
}

bool validateEnSightGoldPattern(const QString& pattern,
                                QString* errorMessage,
                                int lineNumber) {
    if (pattern.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case line %1 contains an empty file name.")
                                .arg(lineNumber);
        }
        return false;
    }
    const QRegularExpression wildcard(QStringLiteral(R"(\*+)"));
    QRegularExpressionMatchIterator matches = wildcard.globalMatch(pattern);
    int groups = 0;
    while (matches.hasNext()) {
        const QRegularExpressionMatch match = matches.next();
        ++groups;
        if (match.capturedLength() > 10) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold case line %1 uses more than ten filename wildcards.")
                                    .arg(lineNumber);
            }
            return false;
        }
    }
    if (groups > 1) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold case line %1 uses multiple disjoint wildcard groups, which the importer cannot expand unambiguously.")
                                .arg(lineNumber);
        }
        return false;
    }
    const QRegularExpression blockWildcard(QStringLiteral(R"(%+)"));
    QRegularExpressionMatchIterator blockMatches = blockWildcard.globalMatch(pattern);
    int blockGroups = 0;
    while (blockMatches.hasNext()) {
        const QRegularExpressionMatch match = blockMatches.next();
        ++blockGroups;
        if (match.capturedLength() > 10) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold case line %1 uses more than ten block-continuation wildcards.")
                                    .arg(lineNumber);
            }
            return false;
        }
    }
    if (blockGroups > 1) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold case line %1 uses multiple disjoint block-continuation wildcard groups.")
                                .arg(lineNumber);
        }
        return false;
    }
    return true;
}

bool appendEnSightGoldReference(EnSightGoldReference reference,
                                EnSightGoldCaseSource* source,
                                QString* errorMessage) {
    if (source == nullptr || reference.patterns.isEmpty()) {
        return false;
    }
    if (source->references.size()
        >= Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "EnSight Gold case exceeds the bounded reference-declaration count.");
        }
        return false;
    }
    for (QString& pattern : reference.patterns) {
        pattern = cleanEnSightGoldName(pattern);
        if (!validateEnSightGoldPattern(pattern, errorMessage, reference.declarationLine)) {
            return false;
        }
    }
    source->references.append(std::move(reference));
    return true;
}

bool addEnSightGoldArray(const QString& name,
                         const QString& declaration,
                         ArrayInfo::Association association,
                          EnSightGoldCaseSource* source,
                          QString* errorMessage,
                          int lineNumber,
                          int componentOverride = 0) {
    if (source == nullptr || name.isEmpty()) {
        return false;
    }
    for (const ArrayInfo& array : std::as_const(source->arrays)) {
        if (array.name.compare(name, Qt::CaseInsensitive) == 0) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold case line %1 redeclares variable name \"%2\" with ambiguous casing or schema.")
                                    .arg(lineNumber)
                                    .arg(name);
            }
            return false;
        }
    }
    ArrayInfo info;
    info.name = name;
    info.components = componentOverride > 0
        ? componentOverride
        : enSightGoldComponentCount(declaration);
    info.association = association;
    source->arrays.append(info);
    return true;
}

bool parseEnSightGoldGeometryDeclaration(const QString& key,
                                         const QString& payload,
                                         int lineNumber,
                                         EnSightGoldCaseSource* source,
                                         QString* errorMessage) {
    QStringList tokens;
    if (!tokenizeEnSightGoldPayload(payload, &tokens, errorMessage, lineNumber)
        || tokens.isEmpty()) {
        if (tokens.isEmpty() && errorMessage != nullptr && errorMessage->isEmpty()) {
            *errorMessage = QStringLiteral("EnSight Gold geometry declaration on line %1 is empty.")
                                .arg(lineNumber);
        }
        return false;
    }
    const QString declaration = key.simplified().toLower();
    if (tokens.contains(QStringLiteral("change_coords_only"), Qt::CaseInsensitive)
        || tokens.contains(QStringLiteral("changing_geometry_per_part"), Qt::CaseInsensitive)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold geometry declaration on line %1 uses a legal changing-geometry variant that this importer does not support safely.")
                                .arg(lineNumber);
        }
        return false;
    }

    EnSightGoldReference reference;
    reference.declaration = declaration;
    reference.declarationLine = lineNumber;
    if (declaration == QStringLiteral("model") || declaration == QStringLiteral("measured")) {
        const bool duplicate = std::any_of(
            source->references.cbegin(), source->references.cend(), [&](const EnSightGoldReference& existing) {
                return existing.role == (declaration == QStringLiteral("model")
                                             ? EnSightGoldReference::Role::ModelGeometry
                                             : EnSightGoldReference::Role::MeasuredGeometry);
            });
        if (duplicate) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold case line %1 duplicates the %2 geometry declaration.")
                                    .arg(lineNumber)
                                    .arg(declaration);
            }
            return false;
        }
        takeEnSightGoldBinding(&tokens, 1, true, &reference.timeSetId, &reference.fileSetId);
        if (tokens.size() != 1) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold %1 geometry on line %2 must contain [time set] [file set] and one quoted-or-unquoted filename.")
                                    .arg(declaration)
                                    .arg(lineNumber);
            }
            return false;
        }
        reference.role = declaration == QStringLiteral("model")
            ? EnSightGoldReference::Role::ModelGeometry
            : EnSightGoldReference::Role::MeasuredGeometry;
        reference.patterns = {tokens.first()};
        return appendEnSightGoldReference(std::move(reference), source, errorMessage);
    }

    static const QSet<QString> auxiliaryDeclarations = {
        QStringLiteral("match"),
        QStringLiteral("boundary"),
        QStringLiteral("rigid_body"),
        QStringLiteral("rigid body"),
        QStringLiteral("vector_glyphs"),
        QStringLiteral("vector glyphs"),
    };
    if (!auxiliaryDeclarations.contains(declaration)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Unsupported EnSight Gold GEOMETRY declaration \"%1\" on line %2.")
                                .arg(key)
                                .arg(lineNumber);
        }
        return false;
    }
    if (tokens.size() != 1) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold auxiliary geometry declaration on line %1 must contain exactly one filename; filenames with spaces must be quoted.")
                                .arg(lineNumber);
        }
        return false;
    }
    reference.role = EnSightGoldReference::Role::Auxiliary;
    reference.patterns = {tokens.first()};
    return appendEnSightGoldReference(std::move(reference), source, errorMessage);
}

bool parseEnSightGoldVariableDeclaration(const QString& key,
                                         const QString& payload,
                                         int lineNumber,
                                         EnSightGoldCaseSource* source,
                                         QString* errorMessage) {
    const QString declaration = key.simplified().toLower();
    QStringList tokens;
    if (!tokenizeEnSightGoldPayload(payload, &tokens, errorMessage, lineNumber)) {
        return false;
    }

    if (declaration == QStringLiteral("constant per case")) {
        int timeSetId = 0;
        if (tokens.size() >= 3) {
            int parsedTimeSet = 0;
            if (parsePositiveIntToken(tokens.first(), &parsedTimeSet)) {
                timeSetId = parsedTimeSet;
                tokens.removeFirst();
            }
        }
        if (tokens.size() < 2) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold constant-per-case declaration on line %1 must contain [time set], a name, and one or more finite numeric values.")
                                    .arg(lineNumber);
            }
            return false;
        }
        QString variableName;
        if (!normalizeEnSightGoldVariableName(tokens.takeFirst(), &variableName,
                                              errorMessage, lineNumber)) {
            return false;
        }
        QVector<double> values;
        QStringList normalizedValues;
        values.reserve(tokens.size());
        normalizedValues.reserve(tokens.size());
        for (const QString& token : std::as_const(tokens)) {
            bool ok = false;
            const double value = token.toDouble(&ok);
            if (!ok || !std::isfinite(value)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Invalid EnSight Gold constant-per-case value on line %1.")
                                        .arg(lineNumber);
                }
                return false;
            }
            values.append(value);
            normalizedValues.append(formatMachineDouble(value));
        }
        if (source->constantsPerCase.size()
            >= Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                    "EnSight Gold case exceeds the bounded constant-per-case declaration count.");
            }
            return false;
        }
        source->constantsPerCase.append(
            EnSightGoldConstantPerCase{lineNumber, timeSetId, variableName, values});
        QStringList normalizedTokens;
        if (timeSetId > 0) {
            normalizedTokens.append(QString::number(timeSetId));
        }
        normalizedTokens.append(variableName);
        normalizedTokens += normalizedValues;
        source->constantPerCaseDeclarations.append(
            QStringLiteral("constant per case: %1")
                .arg(normalizedTokens.join(QLatin1Char(' '))));
        return true;
    }

    EnSightGoldReference reference;
    reference.declaration = declaration;
    reference.declarationLine = lineNumber;
    if (declaration == QStringLiteral("constant per case file")) {
        takeEnSightGoldBinding(&tokens, 2, false, &reference.timeSetId, &reference.fileSetId);
        if (tokens.size() != 2 || tokens.first().isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold constant-per-case-file declaration on line %1 must contain [time set], a name, and one filename.")
                                    .arg(lineNumber);
            }
            return false;
        }
        reference.role = EnSightGoldReference::Role::ConstantVariable;
        if (!normalizeEnSightGoldVariableName(tokens.first(), &reference.variableName,
                                              errorMessage, lineNumber)) {
            return false;
        }
        reference.patterns = {tokens.at(1)};
        return appendEnSightGoldReference(std::move(reference), source, errorMessage);
    }
    if (declaration == QStringLiteral("constant per part")) {
        takeEnSightGoldBinding(&tokens, 2, false, &reference.timeSetId, &reference.fileSetId);
        if ((tokens.size() != 2 && tokens.size() != 3) || tokens.first().isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold constant-per-part declaration on line %1 must contain [time set], a name, a value filename, and an optional index filename.")
                                    .arg(lineNumber);
            }
            return false;
        }
        reference.role = EnSightGoldReference::Role::ConstantVariable;
        if (!normalizeEnSightGoldVariableName(tokens.first(), &reference.variableName,
                                              errorMessage, lineNumber)) {
            return false;
        }
        reference.patterns = tokens.mid(1);
        return appendEnSightGoldReference(std::move(reference), source, errorMessage);
    }

    static const QSet<QString> standardDeclarations = {
        QStringLiteral("scalar per node"),
        QStringLiteral("scalar per element"),
        QStringLiteral("scalar per measured node"),
        QStringLiteral("vector per node"),
        QStringLiteral("vector per element"),
        QStringLiteral("vector per measured node"),
        QStringLiteral("tensor per node"),
        QStringLiteral("tensor per element"),
        QStringLiteral("tensor symm per node"),
        QStringLiteral("tensor symm per element"),
        QStringLiteral("tensor asym per node"),
        QStringLiteral("tensor asym per element"),
    };
    static const QSet<QString> complexDeclarations = {
        QStringLiteral("complex scalar per node"),
        QStringLiteral("complex scalar per element"),
        QStringLiteral("complex vector per node"),
        QStringLiteral("complex vector per element"),
    };
    const bool isStandard = standardDeclarations.contains(declaration);
    const bool isComplex = complexDeclarations.contains(declaration);
    if (!isStandard && !isComplex) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Unsupported EnSight Gold VARIABLE declaration \"%1\" on line %2.")
                                .arg(key)
                                .arg(lineNumber);
        }
        return false;
    }

    takeEnSightGoldBinding(&tokens,
                           isComplex ? 3 : 2,
                           true,
                           &reference.timeSetId,
                           &reference.fileSetId);
    if (isStandard && tokens.size() != 2) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold variable on line %1 must contain [time set] [file set], one name, and one filename; names and paths with spaces must be quoted.")
                                .arg(lineNumber);
        }
        return false;
    }
    if (isComplex && tokens.size() != 4) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold complex variable on line %1 must contain [time set] [file set], one name, real/imaginary filenames, and a finite frequency or UNDEFINED.")
                                .arg(lineNumber);
        }
        return false;
    }
    if (tokens.first().isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold variable name on line %1 is empty.").arg(lineNumber);
        }
        return false;
    }
    if (isComplex && tokens.size() == 4) {
        if (tokens.last().compare(QStringLiteral("UNDEFINED"),
                                  Qt::CaseInsensitive) == 0) {
            reference.frequencyUndefined = true;
            tokens.removeLast();
        } else {
            bool ok = false;
            const double frequency = tokens.last().toDouble(&ok);
            if (!ok || !std::isfinite(frequency)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "EnSight Gold complex frequency on line %1 must be finite numeric data or UNDEFINED.")
                                        .arg(lineNumber);
                }
                return false;
            }
            reference.hasFrequency = true;
            reference.frequency = frequency;
            tokens.removeLast();
        }
    }

    const bool measured = declaration.contains(QStringLiteral("measured"));
    const ArrayInfo::Association association = declaration.contains(QStringLiteral("element"))
        ? ArrayInfo::Association::Cell
        : ArrayInfo::Association::Point;
    QString variableName;
    if (!normalizeEnSightGoldVariableName(tokens.first(), &variableName,
                                          errorMessage, lineNumber)) {
        return false;
    }
    if (!measured) {
        if (isComplex && declaration.contains(QStringLiteral("vector"))) {
            if (!addEnSightGoldArray(variableName + QStringLiteral("_r"), declaration,
                                     association, source, errorMessage,
                                     lineNumber, 3)
                || !addEnSightGoldArray(variableName + QStringLiteral("_i"), declaration,
                                        association, source, errorMessage,
                                        lineNumber, 3)) {
                return false;
            }
        } else if (!addEnSightGoldArray(variableName, declaration, association,
                                        source, errorMessage, lineNumber)) {
            return false;
        }
    }
    reference.role = measured ? EnSightGoldReference::Role::MeasuredVariable
                              : EnSightGoldReference::Role::VolumeVariable;
    reference.variableName = variableName;
    reference.patterns = tokens.mid(1);
    return appendEnSightGoldReference(std::move(reference), source, errorMessage);
}

bool parseEnSightGoldIntegerList(const QString& inlineValues,
                                 int expectedCount,
                                 const QVector<EnSightGoldCaseLine>& lines,
                                 int* lineIndex,
                                 QVector<qint64>* values,
                                 QString* errorMessage,
                                 const QString& label) {
    if (lineIndex == nullptr || values == nullptr || expectedCount <= 0
        || expectedCount > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
        return false;
    }
    values->clear();
    QString current = inlineValues;
    while (true) {
        QStringList tokens;
        if (!tokenizeEnSightGoldPayload(current, &tokens, errorMessage,
                                        lines.at(*lineIndex).number)) {
            return false;
        }
        if (values->size() + tokens.size() > expectedCount) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold %1 contains more than the declared %2 values near line %3.")
                                    .arg(label)
                                    .arg(expectedCount)
                                    .arg(lines.at(*lineIndex).number);
            }
            return false;
        }
        for (const QString& token : tokens) {
            qint64 value = 0;
            if (!parseEnSightGoldIntegerToken(token, &value)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold %1 contains an invalid 32-bit integer on line %2: %3")
                                        .arg(label)
                                        .arg(lines.at(*lineIndex).number)
                                        .arg(token);
                }
                return false;
            }
            values->append(value);
        }
        if (values->size() == expectedCount) {
            return true;
        }
        if (*lineIndex + 1 >= lines.size()) {
            break;
        }
        const EnSightGoldCaseLine& next = lines.at(*lineIndex + 1);
        if (!enSightGoldSectionName(next.text).isEmpty()
            || isEnSightGoldSectionSpelling(next.text)
            || enSightGoldUnquotedColon(next.text) >= 0) {
            break;
        }
        ++(*lineIndex);
        current = next.text;
    }
    if (errorMessage != nullptr) {
        *errorMessage = QStringLiteral("EnSight Gold %1 contains %2 values, but %3 were declared.")
                            .arg(label)
                            .arg(values->size())
                            .arg(expectedCount);
    }
    return false;
}

bool parseEnSightGoldDoubleList(const QString& inlineValues,
                                int expectedCount,
                                const QVector<EnSightGoldCaseLine>& lines,
                                int* lineIndex,
                                QVector<double>* values,
                                QString* errorMessage,
                                const QString& label) {
    if (lineIndex == nullptr || values == nullptr || expectedCount <= 0
        || expectedCount > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
        return false;
    }
    values->clear();
    QString current = inlineValues;
    while (true) {
        QStringList tokens;
        if (!tokenizeEnSightGoldPayload(current, &tokens, errorMessage,
                                        lines.at(*lineIndex).number)) {
            return false;
        }
        if (values->size() + tokens.size() > expectedCount) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold %1 contains more than the declared %2 values near line %3.")
                                    .arg(label)
                                    .arg(expectedCount)
                                    .arg(lines.at(*lineIndex).number);
            }
            return false;
        }
        for (const QString& token : tokens) {
            bool ok = false;
            const double value = token.toDouble(&ok);
            if (!ok || !std::isfinite(value)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold %1 contains a non-finite numeric value on line %2: %3")
                                        .arg(label)
                                        .arg(lines.at(*lineIndex).number)
                                        .arg(token);
                }
                return false;
            }
            values->append(value);
        }
        if (values->size() == expectedCount) {
            return true;
        }
        if (*lineIndex + 1 >= lines.size()) {
            break;
        }
        const EnSightGoldCaseLine& next = lines.at(*lineIndex + 1);
        if (!enSightGoldSectionName(next.text).isEmpty()
            || isEnSightGoldSectionSpelling(next.text)
            || enSightGoldUnquotedColon(next.text) >= 0) {
            break;
        }
        ++(*lineIndex);
        current = next.text;
    }
    if (errorMessage != nullptr) {
        *errorMessage = QStringLiteral("EnSight Gold %1 contains %2 values, but %3 were declared.")
                            .arg(label)
                            .arg(values->size())
                            .arg(expectedCount);
    }
    return false;
}

QString expandEnSightGoldPattern(QString pattern, qint64 filenameNumber) {
    const QRegularExpression stars(QStringLiteral(R"(\*+)"));
    const QRegularExpressionMatch match = stars.match(pattern);
    if (!match.hasMatch()) {
        return pattern;
    }
    const int width = match.capturedLength();
    const QString number = QStringLiteral("%1").arg(filenameNumber, width, 10, QLatin1Char('0'));
    pattern.replace(match.capturedStart(), width, number);
    return pattern;
}

QString expandEnSightGoldBlockPattern(QString pattern, qint64 blockNumber) {
    const QRegularExpression percents(QStringLiteral(R"(%+)"));
    const QRegularExpressionMatch match = percents.match(pattern);
    if (!match.hasMatch()) {
        return pattern;
    }
    const int width = match.capturedLength();
    const QString number = QStringLiteral("%1").arg(blockNumber, width, 10,
                                                    QLatin1Char('0'));
    pattern.replace(match.capturedStart(), width, number);
    return pattern;
}

QString resolveEnSightGoldReferencedPath(const EnSightGoldCaseSource& source,
                                         const QString& expandedPattern) {
    if (QDir::isAbsolutePath(expandedPattern)) {
        return normalizedPath(QFileInfo(expandedPattern).absoluteFilePath());
    }
    return normalizedPath(QFileInfo(QDir(QFileInfo(source.path).absolutePath())
                                        .absoluteFilePath(expandedPattern))
                              .absoluteFilePath());
}

bool finalizeEnSightGoldCaseSource(EnSightGoldCaseSource* source, QString* errorMessage) {
    if (source == nullptr) {
        return false;
    }
    const bool hasAnyBlockContinuation = source->hasBlockContinuationCount
        || source->hasBlockContinuationStart
        || source->hasBlockContinuationIncrement;
    if (hasAnyBlockContinuation
        && (!source->hasBlockContinuationCount
            || !source->hasBlockContinuationStart
            || !source->hasBlockContinuationIncrement
            || source->blockContinuationCount <= 0
            || source->blockContinuationCount
                > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount
            || source->blockContinuationIncrement == 0)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "EnSight Gold BLOCK_CONTINUATION must declare bounded number of sets, filename start number, and non-zero increment exactly once.");
        }
        return false;
    }
    if (hasAnyBlockContinuation) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "EnSight Gold BLOCK_CONTINUATION is legal, but the bundled VTK EnSight reader cannot preserve its structured-block continuation semantics safely; re-export as one merged case, VTU/VTS, or another supported format.");
        }
        return false;
    }

    for (auto timeSetIt = source->timeSets.begin(); timeSetIt != source->timeSets.end(); ++timeSetIt) {
        EnSightGoldTimeSet& timeSet = timeSetIt.value();
        if (!timeSet.hasDeclaredSteps || timeSet.declaredSteps <= 0
            || timeSet.declaredSteps > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount
            || (source->hasMaximumTimeSteps
                && timeSet.declaredSteps > source->maximumTimeSteps)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold time set %1 must declare between 1 and %2 steps.")
                                    .arg(timeSet.id)
                                    .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount);
            }
            return false;
        }
        if (timeSet.timeValues.size() != timeSet.declaredSteps) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold time set %1 has %2 time values, but declares %3 steps.")
                                    .arg(timeSet.id)
                                    .arg(timeSet.timeValues.size())
                                    .arg(timeSet.declaredSteps);
            }
            return false;
        }
        for (qsizetype index = 1; index < timeSet.timeValues.size(); ++index) {
            if (!(timeSet.timeValues.at(index) > timeSet.timeValues.at(index - 1))) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold time set %1 values must be finite and strictly increasing.")
                                        .arg(timeSet.id);
                }
                return false;
            }
        }
        if (timeSet.hasExplicitFilenameNumbers) {
            if (timeSet.hasFilenameStart || timeSet.hasFilenameIncrement
                || timeSet.filenameNumbers.size() != timeSet.declaredSteps) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold time set %1 has conflicting or incomplete filename-number declarations.")
                                        .arg(timeSet.id);
                }
                return false;
            }
        } else if (timeSet.hasFilenameStart || timeSet.hasFilenameIncrement) {
            if (!timeSet.hasFilenameStart || !timeSet.hasFilenameIncrement
                || timeSet.filenameIncrement == 0) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "EnSight Gold time set %1 has an incomplete filename start/increment pair.")
                                        .arg(timeSet.id);
                }
                return false;
            }
            timeSet.filenameNumbers.clear();
            timeSet.filenameNumbers.reserve(timeSet.declaredSteps);
            for (int index = 0; index < timeSet.declaredSteps; ++index) {
                const qint64 number = timeSet.filenameStart
                    + static_cast<qint64>(index) * timeSet.filenameIncrement;
                if (number < std::numeric_limits<int>::min()
                    || number > std::numeric_limits<int>::max()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                                            "EnSight Gold time set %1 filename-number arithmetic exceeds VTK's 32-bit range at step %2.")
                                            .arg(timeSet.id)
                                            .arg(index + 1);
                    }
                    return false;
                }
                timeSet.filenameNumbers.append(number);
            }
        }
    }

    for (const EnSightGoldConstantPerCase& constant : std::as_const(source->constantsPerCase)) {
        if (constant.timeSetId == 0) {
            if (constant.values.size() != 1) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "Static EnSight Gold constant per case on line %1 must contain exactly one value.")
                                        .arg(constant.declarationLine);
                }
                return false;
            }
            continue;
        }
        const auto timeSet = source->timeSets.constFind(constant.timeSetId);
        if (timeSet == source->timeSets.cend()
            || constant.values.size() != timeSet->declaredSteps) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "Transient EnSight Gold constant per case on line %1 must reference an existing time set and provide exactly one value per declared step.")
                                    .arg(constant.declarationLine);
            }
            return false;
        }
    }

    for (auto fileSetIt = source->fileSets.cbegin(); fileSetIt != source->fileSets.cend(); ++fileSetIt) {
        const EnSightGoldFileSet& fileSet = fileSetIt.value();
        qint64 total = 0;
        if (fileSet.chunks.isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold file set %1 does not declare any step chunks.")
                                    .arg(fileSet.id);
            }
            return false;
        }
        for (const EnSightGoldFileChunk& chunk : fileSet.chunks) {
            if (chunk.stepCount <= 0) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold file set %1 contains a non-positive step count.")
                                        .arg(fileSet.id);
                }
                return false;
            }
            total += chunk.stepCount;
            if (total > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold file set %1 exceeds the %2-step limit.")
                                        .arg(fileSet.id)
                                        .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount);
                }
                return false;
            }
        }
    }

    QHash<int, int> fileSetToTimeSet;
    QSet<int> volumeTimeSets;
    source->referencedPaths.clear();
    source->modelGeometryPaths.clear();
    for (EnSightGoldReference& reference : source->references) {
        const EnSightGoldTimeSet* timeSet = nullptr;
        if (reference.timeSetId > 0) {
            const auto timeIt = source->timeSets.constFind(reference.timeSetId);
            if (timeIt == source->timeSets.cend()) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold %1 on line %2 references missing time set %3.")
                                        .arg(reference.declaration)
                                        .arg(reference.declarationLine)
                                        .arg(reference.timeSetId);
                }
                return false;
            }
            timeSet = &timeIt.value();
        }
        const EnSightGoldFileSet* fileSet = nullptr;
        if (reference.fileSetId > 0) {
            if (timeSet == nullptr) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold file set %1 on line %2 is not paired with a time set.")
                                        .arg(reference.fileSetId)
                                        .arg(reference.declarationLine);
                }
                return false;
            }
            const auto fileIt = source->fileSets.constFind(reference.fileSetId);
            if (fileIt == source->fileSets.cend()) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold %1 on line %2 references missing file set %3.")
                                        .arg(reference.declaration)
                                        .arg(reference.declarationLine)
                                        .arg(reference.fileSetId);
                }
                return false;
            }
            fileSet = &fileIt.value();
            const auto boundTime = fileSetToTimeSet.constFind(reference.fileSetId);
            if (boundTime != fileSetToTimeSet.cend() && boundTime.value() != reference.timeSetId) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold file set %1 is bound to multiple time sets (%2 and %3).")
                                        .arg(reference.fileSetId)
                                        .arg(boundTime.value())
                                        .arg(reference.timeSetId);
                }
                return false;
            }
            fileSetToTimeSet.insert(reference.fileSetId, reference.timeSetId);
            const qint64 fileSteps = std::accumulate(
                fileSet->chunks.cbegin(), fileSet->chunks.cend(), qint64{0},
                [](qint64 sum, const EnSightGoldFileChunk& chunk) { return sum + chunk.stepCount; });
            if (fileSteps != timeSet->declaredSteps) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold file set %1 contains %2 steps, but bound time set %3 contains %4.")
                                        .arg(reference.fileSetId)
                                        .arg(fileSteps)
                                        .arg(reference.timeSetId)
                                        .arg(timeSet->declaredSteps);
                }
                return false;
            }
        }
        if ((reference.role == EnSightGoldReference::Role::ModelGeometry
             || reference.role == EnSightGoldReference::Role::VolumeVariable)
            && reference.timeSetId > 0) {
            volumeTimeSets.insert(reference.timeSetId);
        }

        reference.resolvedPaths.clear();
        reference.resolvedPatternPaths.clear();
        reference.expandedPatternNames.clear();
        reference.resolvedPatternPaths.reserve(reference.patterns.size());
        reference.expandedPatternNames.reserve(reference.patterns.size());
        for (const QString& pattern : std::as_const(reference.patterns)) {
            QStringList expandedPatterns;
            const bool hasWildcard = pattern.contains(QLatin1Char('*'));
            const bool hasBlockWildcard = pattern.contains(QLatin1Char('%'));
            const bool participatesInBlockContinuation =
                reference.role != EnSightGoldReference::Role::Auxiliary;
            if (hasBlockWildcard
                && (!hasAnyBlockContinuation || !participatesInBlockContinuation)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                        "EnSight Gold line %1 uses a %% block wildcard outside a supported BLOCK_CONTINUATION geometry/variable declaration.")
                                        .arg(reference.declarationLine);
                }
                return false;
            }
            if (hasAnyBlockContinuation && participatesInBlockContinuation
                && !hasBlockWildcard) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                        "EnSight Gold line %1 must use a %% block wildcard because BLOCK_CONTINUATION is enabled.")
                                        .arg(reference.declarationLine);
                }
                return false;
            }
            if (fileSet != nullptr) {
                if (fileSet->usesFilenameIndices) {
                    if (!hasWildcard) {
                        if (errorMessage != nullptr) {
                            *errorMessage = QStringLiteral("EnSight Gold line %1 binds continuation file set %2 to a filename without wildcards.")
                                                .arg(reference.declarationLine)
                                                .arg(reference.fileSetId);
                        }
                        return false;
                    }
                    for (const EnSightGoldFileChunk& chunk : fileSet->chunks) {
                        expandedPatterns.append(expandEnSightGoldPattern(pattern, chunk.filenameIndex));
                    }
                } else {
                    if (hasWildcard) {
                        if (errorMessage != nullptr) {
                            *errorMessage = QStringLiteral("EnSight Gold line %1 uses wildcards with single-file file set %2, which has no filename indices.")
                                                .arg(reference.declarationLine)
                                                .arg(reference.fileSetId);
                        }
                        return false;
                    }
                    expandedPatterns.append(pattern);
                }
            } else if (timeSet != nullptr && hasWildcard) {
                if (timeSet->filenameNumbers.size() != timeSet->declaredSteps) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "EnSight Gold line %1 uses a transient wildcard, but time set %2 has no complete filename-number declaration.")
                                            .arg(reference.declarationLine)
                                            .arg(reference.timeSetId);
                    }
                    return false;
                }
                for (const qint64 filenameNumber : timeSet->filenameNumbers) {
                    expandedPatterns.append(expandEnSightGoldPattern(pattern, filenameNumber));
                }
            } else {
                if (timeSet == nullptr && hasWildcard) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight Gold line %1 uses a wildcard filename without a time or file set.")
                                            .arg(reference.declarationLine);
                    }
                    return false;
                }
                expandedPatterns.append(pattern);
            }
            if (hasBlockWildcard) {
                QStringList blockExpanded;
                if (expandedPatterns.size()
                    > Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount
                        / source->blockContinuationCount) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "EnSight Gold block/time expansion on line %1 exceeds the bounded dependency count.")
                                            .arg(reference.declarationLine);
                    }
                    return false;
                }
                blockExpanded.reserve(expandedPatterns.size()
                                      * source->blockContinuationCount);
                for (const QString& expanded : std::as_const(expandedPatterns)) {
                    for (int blockIndex = 0;
                         blockIndex < source->blockContinuationCount;
                         ++blockIndex) {
                        const qint64 blockNumber = source->blockContinuationStart
                            + static_cast<qint64>(blockIndex)
                                * source->blockContinuationIncrement;
                        if (blockNumber < std::numeric_limits<int>::min()
                            || blockNumber > std::numeric_limits<int>::max()) {
                            if (errorMessage != nullptr) {
                                *errorMessage = QStringLiteral(
                                    "EnSight Gold BLOCK_CONTINUATION filename arithmetic exceeds the 32-bit range.");
                            }
                            return false;
                        }
                        blockExpanded.append(
                            expandEnSightGoldBlockPattern(expanded, blockNumber));
                    }
                }
                expandedPatterns = blockExpanded;
            }
            QStringList patternPaths;
            patternPaths.reserve(expandedPatterns.size());
            for (const QString& expanded : std::as_const(expandedPatterns)) {
                const QString path = resolveEnSightGoldReferencedPath(*source, expanded);
                patternPaths.append(path);
                if (!containsSourcePath(reference.resolvedPaths, path)) {
                    reference.resolvedPaths.append(path);
                }
                if (!containsSourcePath(source->referencedPaths, path)) {
                    source->referencedPaths.append(path);
                }
                if (reference.role == EnSightGoldReference::Role::ModelGeometry
                    && !containsSourcePath(source->modelGeometryPaths, path)) {
                    source->modelGeometryPaths.append(path);
                }
            }
            reference.resolvedPatternPaths.append(patternPaths);
            reference.expandedPatternNames.append(expandedPatterns);
        }
    }

    if (volumeTimeSets.size() > 1) {
        QStringList ids;
        for (const int id : std::as_const(volumeTimeSets)) {
            ids.append(QString::number(id));
        }
        ids.sort();
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight Gold volume geometry/variables use multiple time sets (%1); Streamcenter+ cannot represent this asynchronous volume timeline safely.")
                                .arg(ids.join(QStringLiteral(", ")));
        }
        return false;
    }
    if (volumeTimeSets.isEmpty()) {
        source->timeValues = {0.0};
    } else {
        source->timeValues = source->timeSets.constFind(*volumeTimeSets.cbegin()).value().timeValues;
    }
    if (source->timeValues.size() > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case exceeds the %1-frame limit.")
                                .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount);
        }
        return false;
    }
    return true;
}

bool parseEnSightGoldCaseSource(const QString& casePath,
                                EnSightGoldCaseSource* source,
                                QString* errorMessage,
                                const QString& logicalCasePath = QString()) {
    if (source == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight Gold metadata storage was not provided.");
        }
        return false;
    }

    const QFileInfo caseInfo(casePath);
    if (!caseInfo.exists() || !caseInfo.isFile() || caseInfo.isSymLink()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("EnSight Gold case file does not exist or is not a regular non-link file: %1")
                                .arg(casePath);
        }
        return false;
    }
    if (caseInfo.suffix().compare(QStringLiteral("case"), Qt::CaseInsensitive) != 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold import expects a *.case file.");
        }
        return false;
    }

    if (caseInfo.size() <= 0
        || caseInfo.size() > Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case metadata must be between 1 and %1 bytes: %2")
                                .arg(Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes)
                                .arg(caseInfo.absoluteFilePath());
        }
        return false;
    }
    QFile file(caseInfo.absoluteFilePath());
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot read EnSight Gold case file: %1").arg(caseInfo.absoluteFilePath());
        }
        return false;
    }

    *source = EnSightGoldCaseSource{};
    source->path = normalizedPath(QFileInfo(logicalCasePath.isEmpty()
                                                ? caseInfo.absoluteFilePath()
                                                : logicalCasePath)
                                      .absoluteFilePath());

    const QByteArray casePayload = file.read(
        Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes + 1);
    if (casePayload.size() > Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes
        || !file.atEnd()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold case metadata exceeds the bounded read limit: %1")
                                .arg(caseInfo.absoluteFilePath());
        }
        return false;
    }
    QVector<EnSightGoldCaseLine> lines;
    if (!readBoundedEnSightGoldCaseLines(casePayload, &lines, errorMessage)) {
        return false;
    }
    QString section;
    bool foundGoldFormat = false;
    QSet<QString> seenSections;
    EnSightGoldTimeSet* currentTimeSet = nullptr;
    EnSightGoldFileSet* currentFileSet = nullptr;
    bool pendingFileSetIndex = false;
    qint64 pendingFilenameIndex = 0;
    for (int lineIndex = 0; lineIndex < lines.size(); ++lineIndex) {
        const EnSightGoldCaseLine& caseLine = lines.at(lineIndex);
        const QString& line = caseLine.text;
        const QString nextSection = enSightGoldSectionName(line);
        if (!nextSection.isEmpty()) {
            if (pendingFileSetIndex) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold FILE section is missing number of steps after filename index near line %1.")
                                        .arg(caseLine.number);
                }
                return false;
            }
            if (nextSection == QStringLiteral("MATERIAL")) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                        "EnSight Gold MATERIAL is a legal optional section, but material/species semantics are not supported safely by this flow-data importer (line %1).")
                                        .arg(caseLine.number);
                }
                return false;
            }
            if (seenSections.contains(nextSection)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold section %1 is declared more than once (line %2).")
                                        .arg(nextSection)
                                        .arg(caseLine.number);
                }
                return false;
            }
            seenSections.insert(nextSection);
            section = nextSection;
            currentTimeSet = nullptr;
            currentFileSet = nullptr;
            continue;
        }
        if (isEnSightGoldSectionSpelling(line)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight Gold section headers must use their required uppercase spelling (line %1): %2")
                                    .arg(caseLine.number)
                                    .arg(line);
            }
            return false;
        }

        if (section == QStringLiteral("QUERY")) {
            // QUERY metadata does not affect geometry, variables, time, or
            // files consumed by VTK.  It is intentionally stripped from the
            // normalized private case after this bounded line parser sees it.
            continue;
        }

        const int colon = enSightGoldUnquotedColon(line);
        if (colon < 0) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Unexpected EnSight Gold content without a declaration key on line %1: %2")
                                    .arg(caseLine.number)
                                    .arg(line);
            }
            return false;
        }

        const QString key = line.left(colon).trimmed();
        const QString normalizedKey = key.simplified().toLower();
        const QString value = line.mid(colon + 1).trimmed();

        if (section == QStringLiteral("FORMAT")) {
            if (normalizedKey != QStringLiteral("type")) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Unsupported EnSight Gold FORMAT declaration \"%1\" on line %2.")
                                        .arg(key)
                                        .arg(caseLine.number);
                }
                return false;
            }
            if (value.simplified().toLower().contains(QStringLiteral("ensight gold"))) {
                foundGoldFormat = true;
            }
            continue;
        }
        if (section == QStringLiteral("GEOMETRY")) {
            if (!parseEnSightGoldGeometryDeclaration(key, value, caseLine.number,
                                                     source, errorMessage)) {
                return false;
            }
            continue;
        }
        if (section == QStringLiteral("VARIABLE")) {
            if (!parseEnSightGoldVariableDeclaration(key, value, caseLine.number,
                                                     source, errorMessage)) {
                return false;
            }
            continue;
        }
        if (section == QStringLiteral("TIME")) {
            if (normalizedKey == QStringLiteral("maximum time steps")) {
                int maximum = 0;
                if (source->hasMaximumTimeSteps || currentTimeSet != nullptr
                    || !source->timeSets.isEmpty()
                    || !parsePositiveIntToken(value, &maximum)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "EnSight Gold maximum time steps on line %1 must be one positive integer declared once before the first time set.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->hasMaximumTimeSteps = true;
                source->maximumTimeSteps = std::min(
                    maximum,
                    static_cast<int>(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount));
                continue;
            }
            if (normalizedKey == QStringLiteral("time set")) {
                QStringList tokens;
                int id = 0;
                if (!tokenizeEnSightGoldPayload(value, &tokens, errorMessage, caseLine.number)
                    || tokens.isEmpty() || !parsePositiveIntToken(tokens.first(), &id)
                    || source->timeSets.contains(id)
                    || source->timeSets.size()
                        >= Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                    if (errorMessage != nullptr && errorMessage->isEmpty()) {
                        *errorMessage = QStringLiteral("Invalid or duplicate EnSight Gold time set on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->timeSets.insert(id, EnSightGoldTimeSet{id, caseLine.number});
                currentTimeSet = &source->timeSets[id];
                continue;
            }
            if (currentTimeSet == nullptr) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight Gold TIME declaration on line %1 appears before a time set.")
                                        .arg(caseLine.number);
                }
                return false;
            }
            if (normalizedKey == QStringLiteral("number of steps")) {
                int steps = 0;
                if (currentTimeSet->hasDeclaredSteps || !parsePositiveIntToken(value, &steps)
                    || steps > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight Gold time set %1 line %2 must declare 1..%3 steps exactly once.")
                                            .arg(currentTimeSet->id)
                                            .arg(caseLine.number)
                                            .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount);
                    }
                    return false;
                }
                currentTimeSet->declaredSteps = steps;
                currentTimeSet->hasDeclaredSteps = true;
            } else if (normalizedKey == QStringLiteral("filename numbers")) {
                if (!currentTimeSet->hasDeclaredSteps || currentTimeSet->hasExplicitFilenameNumbers
                    || currentTimeSet->hasFilenameStart || currentTimeSet->hasFilenameIncrement) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Conflicting or out-of-order filename numbers for EnSight time set %1 on line %2.")
                                            .arg(currentTimeSet->id)
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                currentTimeSet->hasExplicitFilenameNumbers = true;
                if (!parseEnSightGoldIntegerList(value, currentTimeSet->declaredSteps,
                                                 lines, &lineIndex,
                                                 &currentTimeSet->filenameNumbers,
                                                 errorMessage,
                                                 QStringLiteral("filename numbers"))) {
                    return false;
                }
            } else if (normalizedKey == QStringLiteral("filename start number")) {
                qint64 number = 0;
                if (currentTimeSet->hasExplicitFilenameNumbers || currentTimeSet->hasFilenameStart
                    || !parseEnSightGoldIntegerToken(value, &number)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Invalid or duplicate filename start number for EnSight time set %1 on line %2.")
                                            .arg(currentTimeSet->id)
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                currentTimeSet->hasFilenameStart = true;
                currentTimeSet->filenameStart = number;
            } else if (normalizedKey == QStringLiteral("filename increment")) {
                qint64 increment = 0;
                if (currentTimeSet->hasExplicitFilenameNumbers || currentTimeSet->hasFilenameIncrement
                    || !parseEnSightGoldIntegerToken(value, &increment) || increment == 0) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Invalid or duplicate filename increment for EnSight time set %1 on line %2.")
                                            .arg(currentTimeSet->id)
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                currentTimeSet->hasFilenameIncrement = true;
                currentTimeSet->filenameIncrement = increment;
            } else if (normalizedKey == QStringLiteral("time values")) {
                if (!currentTimeSet->hasDeclaredSteps || !currentTimeSet->timeValues.isEmpty()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Invalid, duplicate, or out-of-order time values for EnSight time set %1 on line %2.")
                                            .arg(currentTimeSet->id)
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                if (!parseEnSightGoldDoubleList(value, currentTimeSet->declaredSteps,
                                                lines, &lineIndex,
                                                &currentTimeSet->timeValues,
                                                errorMessage,
                                                QStringLiteral("time values"))) {
                    return false;
                }
            } else if (normalizedKey == QStringLiteral("filename numbers file")
                       || normalizedKey == QStringLiteral("time values file")) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "EnSight Gold external %1 on line %2 is a legal variant that this importer does not support safely.")
                                        .arg(normalizedKey)
                                        .arg(caseLine.number);
                }
                return false;
            } else {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Unsupported EnSight Gold TIME declaration \"%1\" on line %2.")
                                        .arg(key)
                                        .arg(caseLine.number);
                }
                return false;
            }
            continue;
        }
        if (section == QStringLiteral("BLOCK_CONTINUATION")) {
            if (normalizedKey == QStringLiteral("number of sets")) {
                int count = 0;
                if (source->hasBlockContinuationCount
                    || !parsePositiveIntToken(value, &count)
                    || count > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "Invalid or duplicate BLOCK_CONTINUATION number of sets on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->hasBlockContinuationCount = true;
                source->blockContinuationCount = count;
            } else if (normalizedKey == QStringLiteral("filename start number")) {
                qint64 start = 0;
                if (source->hasBlockContinuationStart
                    || !parseEnSightGoldIntegerToken(value, &start)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "Invalid or duplicate BLOCK_CONTINUATION filename start number on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->hasBlockContinuationStart = true;
                source->blockContinuationStart = start;
            } else if (normalizedKey == QStringLiteral("filename increment")) {
                qint64 increment = 0;
                if (source->hasBlockContinuationIncrement
                    || !parseEnSightGoldIntegerToken(value, &increment)
                    || increment == 0) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "Invalid or duplicate BLOCK_CONTINUATION filename increment on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->hasBlockContinuationIncrement = true;
                source->blockContinuationIncrement = increment;
            } else {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                        "Unsupported EnSight Gold BLOCK_CONTINUATION declaration \"%1\" on line %2.")
                                        .arg(key)
                                        .arg(caseLine.number);
                }
                return false;
            }
            continue;
        }
        if (section == QStringLiteral("FILE")) {
            if (normalizedKey == QStringLiteral("file set")) {
                if (pendingFileSetIndex) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight Gold FILE set before line %1 is missing number of steps.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                QStringList tokens;
                int id = 0;
                if (!tokenizeEnSightGoldPayload(value, &tokens, errorMessage, caseLine.number)
                    || tokens.isEmpty() || !parsePositiveIntToken(tokens.first(), &id)
                    || source->fileSets.contains(id)
                    || source->fileSets.size()
                        >= Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                    if (errorMessage != nullptr && errorMessage->isEmpty()) {
                        *errorMessage = QStringLiteral("Invalid or duplicate EnSight Gold file set on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                source->fileSets.insert(id, EnSightGoldFileSet{id, caseLine.number});
                currentFileSet = &source->fileSets[id];
            } else if (normalizedKey == QStringLiteral("filename index")) {
                if (currentFileSet == nullptr || pendingFileSetIndex
                    || (!currentFileSet->chunks.isEmpty() && !currentFileSet->usesFilenameIndices)
                    || !parseEnSightGoldIntegerToken(value, &pendingFilenameIndex)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Invalid filename index in EnSight Gold FILE section on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                currentFileSet->usesFilenameIndices = true;
                pendingFileSetIndex = true;
            } else if (normalizedKey == QStringLiteral("number of steps")) {
                int steps = 0;
                if (currentFileSet == nullptr || !parsePositiveIntToken(value, &steps)
                    || steps > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount
                    || (!pendingFileSetIndex && !currentFileSet->chunks.isEmpty())
                    || (!pendingFileSetIndex && currentFileSet->usesFilenameIndices)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("Invalid or duplicate number of steps in EnSight Gold FILE section on line %1.")
                                            .arg(caseLine.number);
                    }
                    return false;
                }
                currentFileSet->chunks.append(
                    EnSightGoldFileChunk{pendingFileSetIndex, pendingFilenameIndex, steps});
                pendingFileSetIndex = false;
            } else {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Unsupported EnSight Gold FILE declaration \"%1\" on line %2.")
                                        .arg(key)
                                        .arg(caseLine.number);
                }
                return false;
            }
            continue;
        }
        if (section == QStringLiteral("SCRIPTS")) {
            if (normalizedKey != QStringLiteral("metadata")) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Unsupported EnSight Gold SCRIPTS declaration \"%1\" on line %2.")
                                        .arg(key)
                                        .arg(caseLine.number);
                }
                return false;
            }
            QStringList tokens;
            if (!tokenizeEnSightGoldPayload(value, &tokens, errorMessage, caseLine.number)
                || tokens.size() != 1) {
                if (errorMessage != nullptr && errorMessage->isEmpty()) {
                    *errorMessage = QStringLiteral("EnSight Gold metadata script path on line %1 must be one quoted-or-unquoted filename.")
                                        .arg(caseLine.number);
                }
                return false;
            }
            EnSightGoldReference reference;
            reference.role = EnSightGoldReference::Role::Auxiliary;
            reference.declaration = normalizedKey;
            reference.declarationLine = caseLine.number;
            reference.patterns = {tokens.first()};
            if (!appendEnSightGoldReference(std::move(reference), source, errorMessage)) {
                return false;
            }
            continue;
        }
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold declaration appears before a supported section on line %1.")
                                .arg(caseLine.number);
        }
        return false;
    }

    if (pendingFileSetIndex) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold FILE section ends before the final filename index's step count.");
        }
        return false;
    }
    if (!foundGoldFormat) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The selected *.case file is not an EnSight Gold case.");
        }
        return false;
    }
    const bool hasModelGeometry = std::any_of(
        source->references.cbegin(), source->references.cend(), [](const EnSightGoldReference& reference) {
            return reference.role == EnSightGoldReference::Role::ModelGeometry;
        });
    if (!hasModelGeometry) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("The EnSight Gold case does not declare a geometry file: %1").arg(source->path);
        }
        return false;
    }
    if (source->arrays.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("The EnSight Gold case does not declare numeric variables: %1").arg(source->path);
        }
        return false;
    }
    return finalizeEnSightGoldCaseSource(source, errorMessage);
}

bool verifyEnSightGoldReferences(const EnSightGoldCaseSource& source, QString* errorMessage) {
    for (const QString& path : source.referencedPaths) {
        const QFileInfo info(path);
        if (!info.exists() || !info.isFile() || info.isSymLink() || !info.isReadable()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "EnSight Gold case references a missing, unreadable, non-file, or linked dependency: %1")
                                    .arg(path);
            }
            return false;
        }
    }
    return true;
}

bool buildEnSightGoldInventoryFromCases(
    const QVector<EnSightGoldCaseSource>& cases,
    QVector<FlowDataSourceFingerprint>* inventory,
    QString* errorMessage,
    const FlowDataInspector::ProgressCallback& progress,
    const FlowDataInspector::CancellationCallback& shouldCancel) {
    if (inventory == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight inventory output is null.");
        }
        return false;
    }
    inventory->clear();
    QStringList paths;
    QSet<QString> identities;
    qint64 aggregatePathBytes = 0;
    auto appendPath = [&](const QString& path) {
        const QString normalized = normalizedPath(QFileInfo(path).absoluteFilePath());
        const QString identity = sourcePathIdentityKey(normalized);
        if (identities.contains(identity)) {
            return true;
        }
        const qsizetype pathBytes = normalized.toUtf8().size();
        if (paths.size() >= Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount
            || pathBytes > Streamcenter::FlowDataProtocolLimits::MaximumStringBytes
            || aggregatePathBytes + pathBytes
                > Streamcenter::FlowDataProtocolLimits::MaximumAggregateStringBytes) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source inventory exceeds its bounded entry/path limits.");
            }
            return false;
        }
        identities.insert(identity);
        paths.append(normalized);
        aggregatePathBytes += pathBytes;
        return true;
    };

    for (const EnSightGoldCaseSource& source : cases) {
        if (!appendPath(source.path)) {
            return false;
        }
        for (const QString& path : source.referencedPaths) {
            if (!appendPath(path)) {
                return false;
            }
        }
    }
    std::sort(paths.begin(), paths.end(), [](const QString& lhs, const QString& rhs) {
        return sourcePathIdentityKey(lhs) < sourcePathIdentityKey(rhs);
    });

    inventory->reserve(paths.size());
    constexpr qint64 kHashChunkBytes = 4 * 1024 * 1024;
    for (qsizetype index = 0; index < paths.size(); ++index) {
        if (shouldCancel && shouldCancel()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source inventory hashing was canceled.");
            }
            inventory->clear();
            return false;
        }
        const QString& path = paths.at(index);
        const QFileInfo before(path);
        if (!before.exists() || !before.isFile() || before.isSymLink()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source inventory contains a missing, non-file, or linked dependency: %1")
                                    .arg(path);
            }
            inventory->clear();
            return false;
        }
        QFile file(path);
        if (!file.open(QIODevice::ReadOnly)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot hash EnSight source dependency %1: %2")
                                    .arg(path, file.errorString());
            }
            inventory->clear();
            return false;
        }
        QCryptographicHash hash(QCryptographicHash::Sha256);
        qint64 hashedBytes = 0;
        while (!file.atEnd()) {
            if (shouldCancel && shouldCancel()) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("EnSight source inventory hashing was canceled.");
                }
                inventory->clear();
                return false;
            }
            const QByteArray chunk = file.read(kHashChunkBytes);
            if (chunk.isEmpty() && file.error() != QFileDevice::NoError) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral("Cannot finish hashing EnSight dependency %1: %2")
                                        .arg(path, file.errorString());
                }
                inventory->clear();
                return false;
            }
            hash.addData(chunk);
            hashedBytes += chunk.size();
            if (progress) {
                const qint64 totalChunks = std::max<qint64>(1,
                    (before.size() + kHashChunkBytes - 1) / kHashChunkBytes);
                const qint64 completedChunks = std::min(
                    totalChunks,
                    (hashedBytes + kHashChunkBytes - 1) / kHashChunkBytes);
                progress(QStringLiteral("Hashing EnSight source dependency %1 of %2")
                             .arg(index + 1)
                             .arg(paths.size()),
                         static_cast<int>(std::min<qint64>(completedChunks,
                                                          std::numeric_limits<int>::max())),
                         static_cast<int>(std::min<qint64>(totalChunks,
                                                          std::numeric_limits<int>::max())));
            }
        }
        file.close();
        const QFileInfo after(path);
        if (!after.exists() || after.isSymLink()
            || before.size() != after.size()
            || before.lastModified() != after.lastModified()
            || !sameSourcePath(before.canonicalFilePath(), after.canonicalFilePath())) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight dependency changed while it was being hashed: %1").arg(path);
            }
            inventory->clear();
            return false;
        }
        FlowDataSourceFingerprint fingerprint;
        fingerprint.path = sourcePathIdentity(path);
        fingerprint.size = before.size();
        fingerprint.sha256Hex = QString::fromLatin1(hash.result().toHex());
        inventory->append(fingerprint);
        if (progress) {
            progress(QStringLiteral("Hashing EnSight source dependencies"),
                     static_cast<int>(index + 1),
                     static_cast<int>(paths.size()));
        }
    }
    return true;
}

bool enSightGoldSourceLooksLarge(const EnSightGoldCaseSource& source) {
    QStringList resolvedPaths;
    for (const EnSightGoldReference& reference : source.references) {
        for (const QString& path : reference.resolvedPaths) {
            if (!containsSourcePath(resolvedPaths, path)) {
                resolvedPaths.append(path);
            }
        }
    }

    qint64 aggregateBytes = 0;
    for (const QString& path : std::as_const(resolvedPaths)) {
        const qint64 bytes = QFileInfo(path).size();
        if (bytes >= kLargeEnSightGoldSingleFileBytes) {
            return true;
        }
        if (bytes > 0) {
            // Stop at the threshold instead of allowing an untrusted set of
            // file sizes to overflow the signed accumulator.
            if (aggregateBytes >= kLargeEnSightGoldFrameBytes - bytes) {
                return true;
            }
            aggregateBytes += bytes;
        }
    }
    return false;
}

bool anyEnSightGoldSourceLooksLarge(const QVector<EnSightGoldCaseSource>& cases) {
    for (const EnSightGoldCaseSource& source : cases) {
        if (enSightGoldSourceLooksLarge(source)) {
            return true;
        }
    }
    return false;
}

QStringList arrayInfoNames(const QList<ArrayInfo>& arrays) {
    QStringList names;
    for (const ArrayInfo& array : arrays) {
        if (!array.name.trimmed().isEmpty() && !names.contains(array.name, Qt::CaseInsensitive)) {
            names.append(array.name);
        }
    }
    names.sort(Qt::CaseInsensitive);
    return names;
}

bool checkEnSightGoldSourcesForVtk(const QVector<EnSightGoldCaseSource>& cases, QString* errorMessage) {
    Q_UNUSED(cases)
    Q_UNUSED(errorMessage)
    // EnSight Gold binary geometry has several legal encodings (C/Fortran,
    // endian variants, time-step wrappers, ids, structured blocks, nsided,
    // and ghost types).  The former narrow little-endian/no-id parser could
    // misclassify and reject those inputs.  Native parsing is already isolated
    // in the bounded helper, so unknown legal encodings must be delegated to
    // VTK instead of being rejected by an incomplete pre-parser.
    return true;
}

bool timeCodesEquivalent(double lhs, double rhs) {
    return lhs == rhs;
}

QString enSightGoldSourceName(const QVector<EnSightGoldCaseSource>& cases, int caseIndex) {
    if (caseIndex < 0 || caseIndex >= cases.size()) {
        return QStringLiteral("<unknown>");
    }
    return cases.at(caseIndex).path;
}

bool loadEnSightGoldCaseSources(const QStringList& sources,
                                QVector<EnSightGoldCaseSource>* cases,
                                QVector<EnSightGoldFrameSource>* frames,
                                QString* errorMessage) {
    if (cases == nullptr || frames == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight Gold case source storage was not provided.");
        }
        return false;
    }
    cases->clear();
    frames->clear();

    QString collectError;
    const QStringList caseFiles = collectEnSightGoldCaseFiles(sources, &collectError);
    if (caseFiles.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = collectError;
        }
        return false;
    }

    cases->reserve(caseFiles.size());
    qsizetype aggregateFrameCount = 0;
    for (const QString& casePath : caseFiles) {
        EnSightGoldCaseSource source;
        if (!parseEnSightGoldCaseSource(casePath, &source, errorMessage)
            || !verifyEnSightGoldReferences(source, errorMessage)) {
            cases->clear();
            frames->clear();
            return false;
        }
        if (source.timeValues.isEmpty()) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("The EnSight Gold case does not contain any readable frames: %1").arg(casePath);
            }
            cases->clear();
            frames->clear();
            return false;
        }
        if (source.timeValues.size()
                > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount
            || aggregateFrameCount
                    > Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount
                        - source.timeValues.size()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "Selected EnSight Gold cases exceed the %1-frame conversion limit.")
                                    .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount);
            }
            cases->clear();
            frames->clear();
            return false;
        }
        aggregateFrameCount += source.timeValues.size();
        cases->append(source);
    }

    frames->reserve(aggregateFrameCount);
    for (int caseIndex = 0; caseIndex < cases->size(); ++caseIndex) {
        const EnSightGoldCaseSource& source = cases->at(caseIndex);
        for (const double timeValue : source.timeValues) {
            frames->append(EnSightGoldFrameSource{caseIndex, timeValue});
        }
    }

    std::sort(frames->begin(), frames->end(), [](const EnSightGoldFrameSource& lhs, const EnSightGoldFrameSource& rhs) {
        if (lhs.timeValue != rhs.timeValue) {
            return lhs.timeValue < rhs.timeValue;
        }
        return lhs.caseIndex < rhs.caseIndex;
    });
    for (int index = 1; index < frames->size(); ++index) {
        const EnSightGoldFrameSource& previous = frames->at(index - 1);
        const EnSightGoldFrameSource& current = frames->at(index);
        if (!timeCodesEquivalent(previous.timeValue, current.timeValue)) {
            continue;
        }
        if (errorMessage != nullptr) {
            *errorMessage = QString("Duplicate EnSight Gold time code %1 was found in both:\n%2\n%3\nEach selected *.case frame must have a unique time code.")
                                .arg(formatDouble(current.timeValue),
                                     enSightGoldSourceName(*cases, previous.caseIndex),
                                     enSightGoldSourceName(*cases, current.caseIndex));
        }
        return false;
    }
    return true;
}

struct EnSightGoldSnapshotFile {
    FlowDataSourceFingerprint fingerprint;
    QString snapshotPath;
};

bool copyAndVerifySnapshotFile(const QString& sourcePath,
                               const QString& destinationPath,
                               const FlowDataSourceFingerprint& expected,
                               bool requireExpectedSourceIdentity,
                               QString* errorMessage,
                               const FlowDataInspector::CancellationCallback& shouldCancel,
                               const FlowDataInspector::ProgressCallback& progress = {}) {
    const QFileInfo before(sourcePath);
    if (!before.exists() || !before.isFile() || before.isSymLink() || !before.isReadable()
        || before.size() != expected.size
        || (requireExpectedSourceIdentity
            && !sameSourcePath(before.canonicalFilePath(), expected.path))) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight snapshot source is missing, linked, unreadable, has the wrong size, or resolves to a different path: %1")
                                .arg(sourcePath);
        }
        return false;
    }

    QFile source(sourcePath);
    if (!source.open(QIODevice::ReadOnly)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot open EnSight snapshot source %1: %2")
                                .arg(sourcePath, source.errorString());
        }
        return false;
    }
    QSaveFile destination(destinationPath);
    destination.setDirectWriteFallback(false);
    if (!destination.open(QIODevice::WriteOnly)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot create private EnSight snapshot file %1: %2")
                                .arg(destinationPath, destination.errorString());
        }
        return false;
    }

    constexpr qint64 kCopyChunkBytes = 4 * 1024 * 1024;
    QCryptographicHash hash(QCryptographicHash::Sha256);
    qint64 copiedBytes = 0;
    while (!source.atEnd()) {
        if (shouldCancel && shouldCancel()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source snapshot creation was canceled.");
            }
            destination.cancelWriting();
            return false;
        }
        const QByteArray chunk = source.read(kCopyChunkBytes);
        if (chunk.isEmpty() && source.error() != QFileDevice::NoError) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot finish reading EnSight snapshot source %1: %2")
                                    .arg(sourcePath, source.errorString());
            }
            destination.cancelWriting();
            return false;
        }
        if (destination.write(chunk) != chunk.size()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot finish writing private EnSight snapshot file %1: %2")
                                    .arg(destinationPath, destination.errorString());
            }
            destination.cancelWriting();
            return false;
        }
        copiedBytes += chunk.size();
        if (copiedBytes > expected.size) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source grew while it was copied into the snapshot: %1")
                                    .arg(sourcePath);
            }
            destination.cancelWriting();
            return false;
        }
        hash.addData(chunk);
        if (progress) {
            const qint64 totalChunks = std::max<qint64>(
                1, (expected.size + kCopyChunkBytes - 1) / kCopyChunkBytes);
            const qint64 completedChunks = std::min(
                totalChunks,
                (copiedBytes + kCopyChunkBytes - 1) / kCopyChunkBytes);
            progress(QStringLiteral("Copying immutable EnSight dependency"),
                     static_cast<int>(std::min<qint64>(completedChunks,
                                                      std::numeric_limits<int>::max())),
                     static_cast<int>(std::min<qint64>(totalChunks,
                                                      std::numeric_limits<int>::max())));
        }
    }
    source.close();
    const QFileInfo after(sourcePath);
    const QString copiedHash = QString::fromLatin1(hash.result().toHex());
    if (copiedBytes != expected.size || copiedHash != expected.sha256Hex
        || !after.exists() || !after.isFile() || after.isSymLink()
        || after.size() != before.size()
        || after.lastModified() != before.lastModified()
        || !sameSourcePath(after.canonicalFilePath(), before.canonicalFilePath())) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight source dependency changed after validation or did not match the validated inventory while snapshotting: %1")
                                .arg(sourcePath);
        }
        destination.cancelWriting();
        return false;
    }
    if (!destination.commit()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot atomically publish private EnSight snapshot file %1: %2")
                                .arg(destinationPath, destination.errorString());
        }
        return false;
    }
    if (!makeSnapshotFileOwnerReadOnly(destinationPath, errorMessage)) {
        return false;
    }
    const QFileInfo snapshotInfo(destinationPath);
    if (!snapshotInfo.exists() || !snapshotInfo.isFile() || snapshotInfo.isSymLink()
        || snapshotInfo.size() != expected.size) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Private EnSight snapshot file failed its post-copy check: %1")
                                .arg(destinationPath);
        }
        return false;
    }
    return true;
}

bool createVerifiedSnapshotHardLink(const QString& sourcePath,
                                    const QString& destinationPath,
                                    const FlowDataSourceFingerprint& expected,
                                    QString* errorMessage) {
    std::error_code error;
    try {
        std::filesystem::create_hard_link(fsPathFromQString(sourcePath),
                                          fsPathFromQString(destinationPath), error);
    } catch (const std::filesystem::filesystem_error& ex) {
        error = ex.code();
    }
    if (error) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot create private EnSight snapshot hardlink %1: %2")
                                .arg(destinationPath,
                                     QString::fromStdString(error.message()));
        }
        return false;
    }
    bool equivalent = false;
    try {
        equivalent = std::filesystem::equivalent(
            fsPathFromQString(sourcePath), fsPathFromQString(destinationPath), error);
    } catch (const std::filesystem::filesystem_error& ex) {
        error = ex.code();
    }
    const QFileInfo destination(destinationPath);
    if (error || !equivalent || !destination.exists() || !destination.isFile()
        || destination.isSymLink() || isFileSystemReparsePoint(destinationPath)
        || destination.size() != expected.size) {
        QFile::remove(destinationPath);
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "Private EnSight snapshot hardlink failed its identity or safety check: %1")
                                .arg(destinationPath);
        }
        return false;
    }
    return makeSnapshotFileOwnerReadOnly(destinationPath, errorMessage);
}

bool snapshotFileSystemSupportsHardLinks(const QString& verifiedFolder,
                                         const QString& casesFolder) {
    const QString sourcePath = QDir(verifiedFolder).absoluteFilePath(
        QStringLiteral(".hardlink-probe-source"));
    const QString targetPath = QDir(casesFolder).absoluteFilePath(
        QStringLiteral(".hardlink-probe-target"));
    QFile source(sourcePath);
    if (!source.open(QIODevice::WriteOnly) || source.write("probe", 5) != 5) {
        source.close();
        QFile::remove(sourcePath);
        return false;
    }
    source.close();
    std::error_code error;
    bool supported = false;
    try {
        std::filesystem::create_hard_link(fsPathFromQString(sourcePath),
                                          fsPathFromQString(targetPath), error);
        supported = !error
            && std::filesystem::equivalent(fsPathFromQString(sourcePath),
                                           fsPathFromQString(targetPath), error)
            && !error;
    } catch (const std::filesystem::filesystem_error&) {
        supported = false;
    }
    QFile::remove(targetPath);
    QFile::remove(sourcePath);
    return supported;
}

bool addSnapshotRequiredBytes(qint64 value,
                              qint64* total,
                              QString* errorMessage) {
    if (value < 0 || *total > std::numeric_limits<qint64>::max() - value) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight snapshot required-size calculation overflowed.");
        }
        return false;
    }
    *total += value;
    return true;
}

QString snapshotToken(const QString& token) {
    if (token.contains(QRegularExpression(QStringLiteral(R"([\s#,:])")))) {
        return QStringLiteral("\"%1\"").arg(token);
    }
    return token;
}

QString enSightSnapshotPattern(const QString& originalPattern,
                               int referenceIndex,
                               int patternIndex) {
    const QRegularExpression wildcard(QStringLiteral(R"([*%]+)"));
    QString suffix = QFileInfo(originalPattern).suffix().toLower();
    if (!QRegularExpression(QStringLiteral(R"(^[a-z0-9]{1,16}$)"))
             .match(suffix)
             .hasMatch()) {
        suffix = QStringLiteral("bin");
    }
    QString name = QStringLiteral("dependency_%1_%2")
                       .arg(referenceIndex, 5, 10, QLatin1Char('0'))
                       .arg(patternIndex, 2, 10, QLatin1Char('0'));
    QRegularExpressionMatchIterator matches = wildcard.globalMatch(originalPattern);
    while (matches.hasNext()) {
        name += QStringLiteral("_") + matches.next().captured();
    }
    return name + QLatin1Char('.') + suffix;
}

QString remapExpandedSnapshotPattern(const QString& originalPattern,
                                     const QString& mappedPattern,
                                     const QString& expandedPattern,
                                     bool* ok) {
    QString expression = QStringLiteral("^");
    qsizetype cursor = 0;
    QRegularExpressionMatchIterator groups =
        QRegularExpression(QStringLiteral(R"([*%]+)"))
            .globalMatch(originalPattern);
    int groupCount = 0;
    while (groups.hasNext()) {
        const QRegularExpressionMatch group = groups.next();
        expression += QRegularExpression::escape(
            originalPattern.mid(cursor, group.capturedStart() - cursor));
        expression += QStringLiteral("(-?\\d+)");
        cursor = group.capturedEnd();
        ++groupCount;
    }
    expression += QRegularExpression::escape(originalPattern.mid(cursor));
    expression += QLatin1Char('$');
    const QRegularExpressionMatch expansion =
        QRegularExpression(expression).match(expandedPattern);
    if (!expansion.hasMatch()) {
        if (ok != nullptr) {
            *ok = false;
        }
        return {};
    }
    QString result = mappedPattern;
    qsizetype searchOffset = 0;
    for (int groupIndex = 1; groupIndex <= groupCount; ++groupIndex) {
        const QRegularExpressionMatch mappedGroup =
            QRegularExpression(QStringLiteral(R"([*%]+)"))
                .match(result, searchOffset);
        if (!mappedGroup.hasMatch()) {
            if (ok != nullptr) {
                *ok = false;
            }
            return {};
        }
        const QString number = expansion.captured(groupIndex);
        result.replace(mappedGroup.capturedStart(), mappedGroup.capturedLength(), number);
        searchOffset = mappedGroup.capturedStart() + number.size();
    }
    if (ok != nullptr) {
        *ok = true;
    }
    return result;
}

void appendWrappedSnapshotValues(const QString& label,
                                 const QStringList& values,
                                 QStringList* lines) {
    if (lines == nullptr) {
        return;
    }
    lines->append(label + QLatin1Char(':'));
    constexpr qsizetype kMaximumVtkCaseLineBytes = 254;
    QString currentLine;
    for (const QString& value : values) {
        const QString candidate = currentLine.isEmpty()
            ? value
            : currentLine + QLatin1Char(' ') + value;
        if (!currentLine.isEmpty()
            && candidate.toUtf8().size() > kMaximumVtkCaseLineBytes) {
            lines->append(currentLine);
            currentLine = value;
        } else {
            currentLine = candidate;
        }
        if (currentLine.toUtf8().size() > kMaximumVtkCaseLineBytes) {
            // Machine-formatted qint64/double values are far below this, but
            // keep the serializer total if a future token type is added.
            lines->append(currentLine);
            currentLine.clear();
        }
    }
    if (!currentLine.isEmpty()) {
        lines->append(currentLine);
    }
}

bool materializeEnSightGoldSnapshotCase(
    const EnSightGoldCaseSource& source,
    int caseIndex,
    const QString& caseFolder,
    const QHash<QString, EnSightGoldSnapshotFile>& verifiedFiles,
    bool useHardLinks,
    QString* snapshotCasePath,
    QString* errorMessage,
    const FlowDataInspector::ProgressCallback& progress,
    const FlowDataInspector::CancellationCallback& shouldCancel) {
    if (snapshotCasePath == nullptr) {
        return false;
    }
    QVector<QStringList> snapshotPatterns;
    snapshotPatterns.reserve(source.references.size());
    QHash<QString, QString> materializedTargetSources;
    int copiedDependencyCount = 0;
    const int totalDependencyCount = std::accumulate(
        source.references.cbegin(), source.references.cend(), 0,
        [](int sum, const EnSightGoldReference& reference) {
            return sum + std::accumulate(
                             reference.resolvedPatternPaths.cbegin(),
                             reference.resolvedPatternPaths.cend(), 0,
                             [](int patternSum, const QStringList& paths) {
                                 return patternSum + paths.size();
                             });
        });

    for (int referenceIndex = 0; referenceIndex < source.references.size(); ++referenceIndex) {
        const EnSightGoldReference& reference = source.references.at(referenceIndex);
        QStringList mappedPatterns;
        mappedPatterns.reserve(reference.patterns.size());
        for (int patternIndex = 0; patternIndex < reference.patterns.size(); ++patternIndex) {
            const QString& originalPattern = reference.patterns.at(patternIndex);
            const QString mappedPattern = enSightSnapshotPattern(
                originalPattern, referenceIndex, patternIndex);
            mappedPatterns.append(mappedPattern);
            const QStringList& originalPaths = reference.resolvedPatternPaths.at(patternIndex);
            const QStringList& expandedNames =
                reference.expandedPatternNames.at(patternIndex);
            if (expandedNames.size() != originalPaths.size()) {
                if (errorMessage != nullptr) {
                    *errorMessage = QStringLiteral(
                                        "Internal error: EnSight snapshot wildcard expansion count is inconsistent on source line %1.")
                                        .arg(reference.declarationLine);
                }
                return false;
            }
            for (qsizetype pathIndex = 0; pathIndex < originalPaths.size(); ++pathIndex) {
                if (shouldCancel && shouldCancel()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral("EnSight source snapshot creation was canceled.");
                    }
                    return false;
                }
                const QString& originalPath = originalPaths.at(pathIndex);
                const QString identity = sourcePathIdentityKey(originalPath);
                const auto verified = verifiedFiles.constFind(identity);
                if (verified == verifiedFiles.cend()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                                            "Validated EnSight inventory does not contain parsed dependency: %1")
                                            .arg(originalPath);
                    }
                    return false;
                }
                bool mappedOk = false;
                const QString targetName = remapExpandedSnapshotPattern(
                    originalPattern, mappedPattern, expandedNames.at(pathIndex), &mappedOk);
                if (!mappedOk) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                            "Internal error: cannot map expanded EnSight snapshot dependency on source line %1.")
                                            .arg(reference.declarationLine);
                    }
                    return false;
                }
                const QString targetPath = normalizedPath(
                    QDir(caseFolder).absoluteFilePath(targetName));
                const QString targetIdentity = sourcePathIdentityKey(targetPath);
                const auto existingTarget = materializedTargetSources.constFind(
                    targetIdentity);
                if (existingTarget != materializedTargetSources.cend()) {
                    if (existingTarget.value() != identity) {
                        if (errorMessage != nullptr) {
                            *errorMessage = QStringLiteral(
                                "EnSight snapshot maps different dependencies to the same normalized target: %1")
                                                .arg(targetPath);
                        }
                        return false;
                    }
                    ++copiedDependencyCount;
                    if (progress) {
                        progress(QStringLiteral("Reusing immutable EnSight snapshot mapping"),
                                 copiedDependencyCount,
                                 std::max(1, totalDependencyCount));
                    }
                    continue;
                }
                bool materialized = false;
                if (useHardLinks) {
                    QString hardLinkError;
                    materialized = createVerifiedSnapshotHardLink(
                        verified->snapshotPath, targetPath, verified->fingerprint,
                        &hardLinkError);
                    if (!materialized) {
                        constexpr qint64 kFallbackReserve = 16LL * 1024LL * 1024LL;
                        QStorageInfo storage(QFileInfo(targetPath).absolutePath());
                        storage.refresh();
                        if (!storage.isValid() || !storage.isReady()
                            || verified->fingerprint.size
                                > storage.bytesAvailable() - std::min(
                                      storage.bytesAvailable(), kFallbackReserve)) {
                            if (errorMessage != nullptr) {
                                *errorMessage = QStringLiteral(
                                    "EnSight snapshot hardlink failed (%1), and there is insufficient space for the controlled copy fallback of %2 bytes.")
                                                    .arg(hardLinkError)
                                                    .arg(verified->fingerprint.size);
                            }
                            return false;
                        }
                        materialized = copyAndVerifySnapshotFile(
                            verified->snapshotPath, targetPath,
                            verified->fingerprint, false, errorMessage,
                            shouldCancel, progress);
                    }
                } else {
                    materialized = copyAndVerifySnapshotFile(
                        verified->snapshotPath, targetPath,
                        verified->fingerprint, false, errorMessage,
                        shouldCancel, progress);
                }
                if (!materialized) {
                    return false;
                }
                materializedTargetSources.insert(targetIdentity, identity);
                ++copiedDependencyCount;
                if (progress) {
                    progress(QStringLiteral("Materializing immutable EnSight snapshot"),
                             copiedDependencyCount,
                             std::max(1, totalDependencyCount));
                }
            }
        }
        snapshotPatterns.append(mappedPatterns);
    }

    QStringList geometryLines;
    QStringList variableLines;
    for (int constantIndex = 0;
         constantIndex < source.constantsPerCase.size();
         ++constantIndex) {
        if (shouldCancel && shouldCancel()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                    "EnSight source snapshot creation was canceled while materializing constants.");
            }
            return false;
        }
        const EnSightGoldConstantPerCase& constant =
            source.constantsPerCase.at(constantIndex);
        if (constant.timeSetId == 0) {
            variableLines.append(QStringLiteral("constant per case: %1 %2")
                                     .arg(constant.variableName,
                                          formatMachineDouble(
                                              constant.values.constFirst())));
            if (progress) {
                progress(QStringLiteral("Materializing EnSight constants"),
                         constantIndex + 1,
                         source.constantsPerCase.size());
            }
            continue;
        }
        const QString valueFileName = QStringLiteral("constant_%1.values")
                                          .arg(constantIndex, 5, 10,
                                               QLatin1Char('0'));
        const QString valueFilePath = normalizedPath(
            QDir(caseFolder).absoluteFilePath(valueFileName));
        QByteArray valuePayload;
        for (const double value : constant.values) {
            valuePayload += formatMachineDouble(value).toUtf8();
            valuePayload += '\n';
        }
        QSaveFile valueFile(valueFilePath);
        valueFile.setDirectWriteFallback(false);
        if (!valueFile.open(QIODevice::WriteOnly)
            || valueFile.write(valuePayload) != valuePayload.size()
            || !valueFile.commit()
            || !makeSnapshotFileOwnerReadOnly(valueFilePath, errorMessage)) {
            if (errorMessage != nullptr && errorMessage->isEmpty()) {
                *errorMessage = QStringLiteral(
                    "Cannot write transient EnSight constant-per-case snapshot file.");
            }
            return false;
        }
        variableLines.append(
            QStringLiteral("constant per case file: %1 %2 %3")
                .arg(constant.timeSetId)
                .arg(constant.variableName, snapshotToken(valueFileName)));
        if (progress) {
            progress(QStringLiteral("Materializing EnSight transient constants"),
                     constantIndex + 1,
                     source.constantsPerCase.size());
        }
    }
    QStringList scriptLines;
    for (int referenceIndex = 0; referenceIndex < source.references.size(); ++referenceIndex) {
        const EnSightGoldReference& reference = source.references.at(referenceIndex);
        QStringList payload;
        if (reference.timeSetId > 0) {
            payload.append(QString::number(reference.timeSetId));
        }
        if (reference.fileSetId > 0) {
            payload.append(QString::number(reference.fileSetId));
        }
        const QStringList& mappedPatterns = snapshotPatterns.at(referenceIndex);
        if (reference.role == EnSightGoldReference::Role::ModelGeometry
            || reference.role == EnSightGoldReference::Role::MeasuredGeometry) {
            payload += mappedPatterns;
            geometryLines.append(reference.declaration + QStringLiteral(": ")
                                 + payload.join(QLatin1Char(' ')));
        } else if (reference.role == EnSightGoldReference::Role::Auxiliary) {
            payload += mappedPatterns;
            const QString line = reference.declaration + QStringLiteral(": ")
                + payload.join(QLatin1Char(' '));
            if (reference.declaration == QStringLiteral("metadata")) {
                scriptLines.append(line);
            } else {
                geometryLines.append(line);
            }
        } else {
            payload.append(snapshotToken(reference.variableName));
            payload += mappedPatterns;
            if (reference.hasFrequency) {
                payload.append(formatMachineDouble(reference.frequency));
            } else if (reference.frequencyUndefined) {
                payload.append(QStringLiteral("UNDEFINED"));
            }
            variableLines.append(reference.declaration + QStringLiteral(": ")
                                 + payload.join(QLatin1Char(' ')));
        }
    }

    QStringList outputLines = {
        QStringLiteral("FORMAT"),
        QStringLiteral("type: ensight gold"),
        QString(),
        QStringLiteral("GEOMETRY"),
    };
    outputLines += geometryLines;
    outputLines.append(QString());
    outputLines.append(QStringLiteral("VARIABLE"));
    outputLines += variableLines;

    if (!source.timeSets.isEmpty()) {
        outputLines.append(QString());
        outputLines.append(QStringLiteral("TIME"));
        if (source.hasMaximumTimeSteps) {
            outputLines.append(QStringLiteral("maximum time steps: %1")
                                   .arg(source.maximumTimeSteps));
        }
        QList<int> ids = source.timeSets.keys();
        std::sort(ids.begin(), ids.end());
        for (const int id : std::as_const(ids)) {
            const EnSightGoldTimeSet& timeSet = source.timeSets[id];
            outputLines.append(QStringLiteral("time set: %1").arg(id));
            outputLines.append(QStringLiteral("number of steps: %1")
                                   .arg(timeSet.declaredSteps));
            QStringList filenameNumbers;
            filenameNumbers.reserve(timeSet.filenameNumbers.size());
            for (const qint64 number : timeSet.filenameNumbers) {
                filenameNumbers.append(QString::number(number));
            }
            if (!filenameNumbers.isEmpty()) {
                appendWrappedSnapshotValues(QStringLiteral("filename numbers"),
                                            filenameNumbers, &outputLines);
            }
            QStringList timeValues;
            timeValues.reserve(timeSet.timeValues.size());
            for (const double value : timeSet.timeValues) {
                timeValues.append(formatMachineDouble(value));
            }
            appendWrappedSnapshotValues(QStringLiteral("time values"),
                                        timeValues, &outputLines);
        }
    }
    if (!source.fileSets.isEmpty()) {
        outputLines.append(QString());
        outputLines.append(QStringLiteral("FILE"));
        QList<int> ids = source.fileSets.keys();
        std::sort(ids.begin(), ids.end());
        for (const int id : std::as_const(ids)) {
            const EnSightGoldFileSet& fileSet = source.fileSets[id];
            outputLines.append(QStringLiteral("file set: %1").arg(id));
            for (const EnSightGoldFileChunk& chunk : fileSet.chunks) {
                if (chunk.hasFilenameIndex) {
                    outputLines.append(QStringLiteral("filename index: %1")
                                           .arg(chunk.filenameIndex));
                }
                outputLines.append(QStringLiteral("number of steps: %1")
                                       .arg(chunk.stepCount));
            }
        }
    }
    if (!scriptLines.isEmpty()) {
        outputLines.append(QString());
        outputLines.append(QStringLiteral("SCRIPTS"));
        outputLines += scriptLines;
    }
    if (source.hasBlockContinuationCount) {
        outputLines.append(QString());
        outputLines.append(QStringLiteral("BLOCK_CONTINUATION"));
        outputLines.append(QStringLiteral("number of sets: %1")
                               .arg(source.blockContinuationCount));
        outputLines.append(QStringLiteral("filename start number: %1")
                               .arg(source.blockContinuationStart));
        outputLines.append(QStringLiteral("filename increment: %1")
                               .arg(source.blockContinuationIncrement));
    }

    for (const QString& line : std::as_const(outputLines)) {
        if (line.toUtf8().size() > 254) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                                    "Normalized EnSight snapshot case would exceed the bundled VTK reader's 254-byte safe line limit.");
            }
            return false;
        }
    }
    const QByteArray payload = outputLines.join(QLatin1Char('\n')).toUtf8()
        + QByteArrayLiteral("\n");
    if (payload.size() > Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Normalized EnSight snapshot case exceeds the bounded metadata size.");
        }
        return false;
    }
    const QString casePath = normalizedPath(
        QDir(caseFolder).absoluteFilePath(
            QStringLiteral("snapshot_%1.case").arg(caseIndex, 5, 10, QLatin1Char('0'))));
    QSaveFile caseFile(casePath);
    caseFile.setDirectWriteFallback(false);
    if (!caseFile.open(QIODevice::WriteOnly)
        || caseFile.write(payload) != payload.size()
        || !caseFile.commit()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot write normalized private EnSight snapshot case: %1")
                                .arg(casePath);
        }
        return false;
    }
    if (!makeSnapshotFileOwnerReadOnly(casePath, errorMessage)) {
        return false;
    }
    *snapshotCasePath = casePath;
    return true;
}

}  // namespace

bool FlowDataInspector::buildEnSightGoldSourceInventory(
    const QStringList& casePaths,
    QVector<FlowDataSourceFingerprint>* inventory,
    QString* errorMessage,
    const ProgressCallback& progress,
    const CancellationCallback& shouldCancel) {
    QVector<EnSightGoldCaseSource> cases;
    QVector<EnSightGoldFrameSource> frames;
    if (!loadEnSightGoldCaseSources(casePaths, &cases, &frames, errorMessage)) {
        return false;
    }
    return buildEnSightGoldInventoryFromCases(
        cases, inventory, errorMessage, progress, shouldCancel);
}

bool FlowDataInspector::materializeEnSightGoldSourceSnapshot(
    const QStringList& casePaths,
    const QVector<FlowDataSourceFingerprint>& expectedInventory,
    const QString& snapshotRoot,
    QStringList* snapshotCasePaths,
    QString* errorMessage,
    const ProgressCallback& progress,
    const CancellationCallback& shouldCancel) {
    if (snapshotCasePaths == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight snapshot case output is null.");
        }
        return false;
    }
    snapshotCasePaths->clear();
    const QFileInfo rootInfo(snapshotRoot);
    if (snapshotRoot.trimmed().isEmpty() || !QDir::isAbsolutePath(snapshotRoot)
        || !rootInfo.exists() || !rootInfo.isDir() || rootInfo.isSymLink()
        || directoryChainContainsReparsePoint(snapshotRoot)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "EnSight snapshot root must be an existing absolute regular, non-link directory.");
        }
        return false;
    }
    QDir root(snapshotRoot);
    if (!root.entryList(QDir::AllEntries | QDir::NoDotAndDotDot
                        | QDir::Hidden | QDir::System)
             .isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight snapshot root must be private and empty before materialization.");
        }
        return false;
    }
    if (expectedInventory.isEmpty()
        || expectedInventory.size()
            > Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Validated EnSight source inventory is empty or exceeds its entry limit.");
        }
        return false;
    }
#ifdef Q_OS_UNIX
    if (!QFile::setPermissions(snapshotRoot,
                               QFileDevice::ReadOwner | QFileDevice::WriteOwner
                                   | QFileDevice::ExeOwner)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot restrict EnSight snapshot root permissions to its owner.");
        }
        return false;
    }
#endif

    const QString verifiedFolder = normalizedPath(
        root.absoluteFilePath(QStringLiteral("verified")));
    const QString casesFolder = normalizedPath(
        root.absoluteFilePath(QStringLiteral("cases")));
    if (!QDir().mkpath(verifiedFolder) || !QDir().mkpath(casesFolder)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot create private EnSight snapshot subdirectories.");
        }
        return false;
    }
    if (directoryChainContainsReparsePoint(verifiedFolder)
        || directoryChainContainsReparsePoint(casesFolder)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "EnSight snapshot subdirectories must not resolve through links or reparse points.");
        }
        return false;
    }
#ifdef Q_OS_UNIX
    const auto ownerDirectoryPermissions = QFileDevice::ReadOwner
        | QFileDevice::WriteOwner | QFileDevice::ExeOwner;
    if (!QFile::setPermissions(verifiedFolder, ownerDirectoryPermissions)
        || !QFile::setPermissions(casesFolder, ownerDirectoryPermissions)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot restrict EnSight snapshot subdirectory permissions.");
        }
        return false;
    }
#endif

    const bool useHardLinks = snapshotFileSystemSupportsHardLinks(
        verifiedFolder, casesFolder);
    constexpr qint64 kSnapshotFreeSpaceReserve = 16LL * 1024LL * 1024LL;
    qint64 initialRequiredBytes = kSnapshotFreeSpaceReserve;
    for (const FlowDataSourceFingerprint& fingerprint : expectedInventory) {
        if (!addSnapshotRequiredBytes(fingerprint.size,
                                      &initialRequiredBytes, errorMessage)) {
            return false;
        }
    }
    QStorageInfo snapshotStorage(snapshotRoot);
    snapshotStorage.refresh();
    if (!snapshotStorage.isValid() || !snapshotStorage.isReady()
        || snapshotStorage.bytesAvailable() < initialRequiredBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "Insufficient free space for the verified EnSight snapshot: %1 bytes are required, %2 bytes are available.")
                                .arg(initialRequiredBytes)
                                .arg(snapshotStorage.bytesAvailable());
        }
        return false;
    }

    QHash<QString, EnSightGoldSnapshotFile> verifiedFiles;
    verifiedFiles.reserve(expectedInventory.size());
    for (qsizetype index = 0; index < expectedInventory.size(); ++index) {
        if (shouldCancel && shouldCancel()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("EnSight source snapshot creation was canceled.");
            }
            return false;
        }
        const FlowDataSourceFingerprint& fingerprint = expectedInventory.at(index);
        const QString identity = sourcePathIdentityKey(fingerprint.path);
        if (fingerprint.path.isEmpty() || !QDir::isAbsolutePath(fingerprint.path)
            || fingerprint.size < 0 || fingerprint.sha256Hex.size() != 64
            || verifiedFiles.contains(identity)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Validated EnSight inventory contains an invalid or duplicate fingerprint.");
            }
            return false;
        }
        const QString suffix = QFileInfo(fingerprint.path).suffix()
                                   .compare(QStringLiteral("case"), Qt::CaseInsensitive) == 0
            ? QStringLiteral("case")
            : QStringLiteral("src");
        const QString copyPath = normalizedPath(
            QDir(verifiedFolder).absoluteFilePath(
                QStringLiteral("source_%1.%2")
                    .arg(index, 5, 10, QLatin1Char('0'))
                    .arg(suffix)));
        if (!copyAndVerifySnapshotFile(fingerprint.path,
                                       copyPath,
                                       fingerprint,
                                       true,
                                       errorMessage,
                                       shouldCancel,
                                       progress)) {
            return false;
        }
        verifiedFiles.insert(identity, EnSightGoldSnapshotFile{fingerprint, copyPath});
        if (progress) {
            progress(QStringLiteral("Verifying EnSight source snapshot"),
                     static_cast<int>(index + 1),
                     static_cast<int>(expectedInventory.size()));
        }
    }

    QString collectError;
    const QStringList actualCasePaths = collectEnSightGoldCaseFiles(casePaths, &collectError);
    if (actualCasePaths.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = collectError;
        }
        return false;
    }
    QVector<EnSightGoldCaseSource> parsedCases;
    parsedCases.reserve(actualCasePaths.size());
    QSet<QString> parsedInventoryIdentities;
    qsizetype aggregateReferenceCount = 0;
    for (const QString& requestedCasePath : actualCasePaths) {
        const QString logicalCasePath = sourcePathIdentity(requestedCasePath);
        const QString identity = sourcePathIdentityKey(logicalCasePath);
        const auto copiedCase = verifiedFiles.constFind(identity);
        if (copiedCase == verifiedFiles.cend()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Validated EnSight inventory does not contain selected case: %1")
                                    .arg(requestedCasePath);
            }
            return false;
        }
        EnSightGoldCaseSource parsed;
        if (!parseEnSightGoldCaseSource(copiedCase->snapshotPath,
                                        &parsed,
                                        errorMessage,
                                        logicalCasePath)) {
            return false;
        }
        parsedInventoryIdentities.insert(sourcePathIdentityKey(parsed.path));
        for (const QString& dependency : parsed.referencedPaths) {
            parsedInventoryIdentities.insert(sourcePathIdentityKey(dependency));
        }
        for (const EnSightGoldReference& reference : parsed.references) {
            for (const QStringList& paths : reference.resolvedPatternPaths) {
                if (aggregateReferenceCount
                    > Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount
                        - paths.size()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QStringLiteral(
                                            "EnSight snapshot expansion exceeds the bounded dependency-reference count.");
                    }
                    return false;
                }
                aggregateReferenceCount += paths.size();
            }
        }
        if (aggregateReferenceCount
            > Streamcenter::FlowDataProtocolLimits::MaximumInventoryEntryCount
                - parsed.constantsPerCase.size()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral(
                    "EnSight snapshot exceeds the unified dependency/generated-file entry budget.");
            }
            return false;
        }
        aggregateReferenceCount += parsed.constantsPerCase.size();
        parsedCases.append(std::move(parsed));
    }
    QSet<QString> expectedIdentities;
    for (const FlowDataSourceFingerprint& fingerprint : expectedInventory) {
        expectedIdentities.insert(sourcePathIdentityKey(fingerprint.path));
    }
    if (parsedInventoryIdentities != expectedIdentities) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "Validated EnSight inventory does not exactly match dependencies parsed from the immutable case snapshot.");
        }
        return false;
    }

    qint64 remainingRequiredBytes = kSnapshotFreeSpaceReserve;
    if (!addSnapshotRequiredBytes(
            static_cast<qint64>(parsedCases.size())
                * Streamcenter::FlowDataProtocolLimits::MaximumCaseMetadataBytes,
            &remainingRequiredBytes, errorMessage)) {
        return false;
    }
    if (!useHardLinks) {
        for (const EnSightGoldCaseSource& parsed : std::as_const(parsedCases)) {
            for (const EnSightGoldReference& reference : parsed.references) {
                for (const QStringList& paths : reference.resolvedPatternPaths) {
                    for (const QString& path : paths) {
                        const auto verified = verifiedFiles.constFind(
                            sourcePathIdentityKey(path));
                        if (verified == verifiedFiles.cend()) {
                            if (errorMessage != nullptr) {
                                *errorMessage = QStringLiteral(
                                    "Validated EnSight inventory is missing a dependency during space preflight: %1")
                                                    .arg(path);
                            }
                            return false;
                        }
                        if (!addSnapshotRequiredBytes(verified->fingerprint.size,
                                                      &remainingRequiredBytes,
                                                      errorMessage)) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    snapshotStorage.refresh();
    if (snapshotStorage.bytesAvailable() < remainingRequiredBytes) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                                "Insufficient free space to materialize EnSight snapshot mappings: %1 bytes are required, %2 bytes are available.")
                                .arg(remainingRequiredBytes)
                                .arg(snapshotStorage.bytesAvailable());
        }
        return false;
    }

    snapshotCasePaths->reserve(parsedCases.size());
    for (int caseIndex = 0; caseIndex < parsedCases.size(); ++caseIndex) {
        const QString caseFolder = normalizedPath(
            QDir(casesFolder).absoluteFilePath(
                QStringLiteral("case_%1").arg(caseIndex, 5, 10, QLatin1Char('0'))));
        if (!QDir().mkpath(caseFolder)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot create private EnSight snapshot case directory.");
            }
            snapshotCasePaths->clear();
            return false;
        }
#ifdef Q_OS_UNIX
        if (!QFile::setPermissions(caseFolder,
                                   QFileDevice::ReadOwner | QFileDevice::WriteOwner
                                       | QFileDevice::ExeOwner)) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot restrict EnSight snapshot case directory permissions.");
            }
            snapshotCasePaths->clear();
            return false;
        }
#endif
        QString snapshotCase;
        if (!materializeEnSightGoldSnapshotCase(parsedCases.at(caseIndex),
                                                caseIndex,
                                                 caseFolder,
                                                 verifiedFiles,
                                                 useHardLinks,
                                                 &snapshotCase,
                                                errorMessage,
                                                progress,
                                                shouldCancel)) {
            snapshotCasePaths->clear();
            return false;
        }
        snapshotCasePaths->append(snapshotCase);
    }
    return !snapshotCasePaths->isEmpty();
}

FlowDataScanResult FlowDataInspector::scanFolder(const QString& folder,
                                                 FlowDataKind dataKind,
                                                 bool coordinateCheck,
                                                 const ProgressCallback& progress,
                                                 const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    const QFileInfo sourceInfo(folder);
    const bool singleFile = sourceInfo.exists() && sourceInfo.isFile();
    if (singleFile && isNek5000DescriptorExtension(sourceInfo.suffix())) {
        // Keep the generic FlowData entry point on the same canonical native
        // reader path as the dedicated Nek5000 import flow.
        return scanNek5000File(sourceInfo.absoluteFilePath(), dataKind,
                               coordinateCheck, progress, shouldCancel);
    }
    result.sourceFolder = normalizedPath(singleFile ? sourceInfo.absoluteFilePath() : folder);
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataName = singleFile ? sourceInfo.completeBaseName() : QFileInfo(folder).fileName();
    result.entry.dataVersion = QDate::currentDate().toString(Qt::ISODate);
    result.entry.fileType = QStringLiteral("vtk series");

    auto canceled = [&]() {
        if (shouldCancel && shouldCancel()) {
            result.errorMessage = QStringLiteral("Flow data validation interrupted.");
            return true;
        }
        return false;
    };
    if (canceled()) {
        return result;
    }

    if (!singleFile) {
        const QDir dir(folder);
        if (!dir.exists()) {
            result.errorMessage = QString("Flow data folder does not exist: %1").arg(folder);
            return result;
        }
    } else if (!sourceInfo.isReadable()) {
        result.errorMessage = QString("Flow data file is not readable: %1").arg(folder);
        return result;
    }

    result.files = collectVtkFiles(folder, dataKind == FlowDataKind::Raw);
    if (result.files.isEmpty()) {
        const QString extensions = dataKind == FlowDataKind::Raw
            ? QStringLiteral("*.vtk, *.vti, *.vtr, *.vts, *.vtu, *.vtm, *.vtkhdf, *.dat, *.plt, *.nek5000, *.nek")
            : QStringLiteral("*.vtk, *.vti, *.vtr, *.vts, *.vtu, *.vtm, *.vtkhdf, *.dat, *.plt");
        result.errorMessage = QString("No supported data files (%1) were found in %2.").arg(extensions, folder);
        return result;
    }
    const bool onlyNek5000Descriptors = std::all_of(
        result.files.cbegin(), result.files.cend(), [](const QString& file) {
            return isNek5000DescriptorExtension(vtkExtensionSuffix(file));
        });
    const bool containsNek5000Descriptor = std::any_of(
        result.files.cbegin(), result.files.cend(), [](const QString& file) {
            return isNek5000DescriptorExtension(vtkExtensionSuffix(file));
        });
    if (onlyNek5000Descriptors) {
        // .nek and .nek5000 are aliases for the same descriptor.  A folder
        // with a single native series should behave exactly like direct
        // import, while several independent descriptors require the caller to
        // choose one rather than silently treating descriptors as time steps.
        QStringList descriptors;
        const auto sameDescriptorSeries = [](const QFileInfo& lhs,
                                             const QFileInfo& rhs) {
            return sameSourcePath(lhs.absolutePath(), rhs.absolutePath())
                && lhs.completeBaseName().compare(rhs.completeBaseName(),
                                                   fileSystemPathCaseSensitivity()) == 0;
        };
        for (const QString& file : std::as_const(result.files)) {
            const QFileInfo info(file);
            int existingIndex = -1;
            for (int index = 0; index < descriptors.size(); ++index) {
                const QFileInfo existing(descriptors.at(index));
                if (sameDescriptorSeries(existing, info)) {
                    existingIndex = index;
                    break;
                }
            }
            if (existingIndex < 0) {
                descriptors.append(info.absoluteFilePath());
            } else if (info.suffix().compare(QStringLiteral("nek5000"),
                                              Qt::CaseInsensitive) == 0) {
                // Prefer the canonical extension when both aliases exist.
                descriptors[existingIndex] = info.absoluteFilePath();
            }
        }
        if (descriptors.size() != 1) {
            result.errorMessage = QStringLiteral(
                "A FlowData folder contains multiple Nek5000 descriptors. "
                "Select one *.nek5000 or *.nek descriptor explicitly.");
            return result;
        }
        return scanNek5000File(descriptors.front(), dataKind,
                               coordinateCheck, progress, shouldCancel);
    }
    if (containsNek5000Descriptor) {
        result.errorMessage = QStringLiteral(
            "A FlowData source cannot mix a Native Nek5000 descriptor with generic VTK files. "
            "Select the Nek5000 descriptor explicitly.");
        return result;
    }
    if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
            return vtkExtensionSuffix(file).compare(QStringLiteral("dat"), Qt::CaseInsensitive) == 0;
        })) {
        result.entry.fileType = QStringLiteral("tecplot ASCII series");
    } else if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
                   return vtkExtensionSuffix(file).compare(QStringLiteral("plt"), Qt::CaseInsensitive) == 0;
               })) {
        result.entry.fileType = QStringLiteral("tecplot binary series");
    } else if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
                   return isTecplotExtension(vtkExtensionSuffix(file));
               })) {
        result.entry.fileType = QStringLiteral("tecplot series");
    } else if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
            const QString suffix = vtkExtensionSuffix(file);
            return suffix.compare(QStringLiteral("nek5000"), Qt::CaseInsensitive) == 0 ||
                   suffix.compare(QStringLiteral("nek"), Qt::CaseInsensitive) == 0;
        })) {
        result.entry.fileType = QStringLiteral("nek5000 series");
    } else if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
                   return vtkExtensionSuffix(file).compare(QStringLiteral("vtkhdf"), Qt::CaseInsensitive) == 0;
               })) {
        result.entry.fileType = QStringLiteral("vtkhdf series");
    }
    const bool validateEveryFrame = coordinateCheck
        || dataKind == FlowDataKind::Resampled;
    if (progress) {
        progress(QStringLiteral("Validating flow data"), 1,
                 validateEveryFrame ? result.files.size() : 2);
    }
    if (canceled()) {
        return result;
    }

    QString readError;
    VtkDataSummary first;
    if (!summarizeAnyVtk(result.files.first(), &first, &readError)) {
        result.errorMessage = readError;
        return result;
    }

    int dims[3] = {0, 0, 0};
    if (first.hasDims) {
        for (int index = 0; index < 3; ++index) {
            dims[index] = first.dims[index];
        }
    }
    if (dataKind == FlowDataKind::Resampled
        && (!first.uniformStructured || !first.hasDims)) {
        result.errorMessage = "Resampled flow data must be a uniform structured VTK dataset.";
        return result;
    }

    const double* bounds = first.bounds;
    result.entry.xRange = formatRange(bounds[0], bounds[1]);
    result.entry.yRange = formatRange(bounds[2], bounds[3]);
    result.entry.zRange = formatRange(bounds[4], bounds[5]);
    const int dimension = first.spatialDimension;
    result.entry.dimension = QString::number(dimension);

    if (dataKind == FlowDataKind::Resampled) {
        result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                     .arg(std::max(1, dims[0]))
                                     .arg(std::max(1, dims[1]))
                                     .arg(std::max(1, dims[2]));
    } else {
        result.entry.meshCount = QString::number(result.files.size());
    }

    result.termCandidates = first.termCandidates;
    inferTerms(first.arrays, dimension, &result.entry);
    result.entry.termCandidates = result.termCandidates;

    if (result.termCandidates.isEmpty()) {
        result.errorMessage = "The first data snapshot does not contain numeric point or cell arrays.";
        return result;
    }

    if (validateEveryFrame) {
        for (int index = 1; index < result.files.size(); ++index) {
            if (canceled()) {
                return result;
            }
            QString frameError;
            VtkDataSummary frame;
            const bool frameOk = summarizeAnyVtk(result.files.at(index), &frame, &frameError);
            if (canceled()) {
                return result;
            }
            if (!frameOk) {
                result.errorMessage = frameError;
                return result;
            }
            if (dataKind == FlowDataKind::Resampled
                && (!frame.uniformStructured || !frame.hasDims)) {
                result.errorMessage = QString(
                    "Every resampled flow data frame must be a uniform structured VTK dataset: %1")
                                          .arg(result.files.at(index));
                return result;
            }
            if (dataKind == FlowDataKind::Resampled
                && (first.dims[0] != frame.dims[0]
                    || first.dims[1] != frame.dims[1]
                    || first.dims[2] != frame.dims[2])) {
                result.errorMessage = QString(
                    "Resampled flow data grid dimensions differ from the first frame: %1")
                                          .arg(result.files.at(index));
                return result;
            }
            if (dataKind == FlowDataKind::Resampled
                && first.spatialDimension != frame.spatialDimension) {
                result.errorMessage = QString(
                    "Resampled flow data spatial dimension differs from the first frame: %1")
                                          .arg(result.files.at(index));
                return result;
            }
            if (dataKind == FlowDataKind::Resampled
                && (!sameCandidateFieldSet(first.termCandidates,
                                           frame.termCandidates)
                    || !sameArraySchema(first.arrays, frame.arrays, true))) {
                result.errorMessage = QString(
                    "Resampled flow data field/component schema differs from the first frame: %1")
                                          .arg(result.files.at(index));
                return result;
            }
            if (coordinateCheck && !coordinatesEquivalent(first, frame)) {
                result.errorMessage = QString("Snapshot grid coordinates differ from the first frame: %1")
                                          .arg(result.files.at(index));
                return result;
            }
            if (progress) {
                progress(QStringLiteral("Validating flow data"), index + 1, result.files.size());
            }
            if (canceled()) {
                return result;
            }
        }
    } else if (progress) {
        progress(QStringLiteral("Validating flow data"), 2, 2);
        if (canceled()) {
            return result;
        }
    }

    result.entry.checkStatus = QStringLiteral("passed");
    result.passed = true;
    return result;
}

bool fingerprintStableRegularFile(
    const QString& filePath,
    const QFileInfo& expectedInfo,
    FlowDataSourceFingerprint* fingerprint,
    QString* errorMessage,
    const FlowDataInspector::ProgressCallback& progress,
    const FlowDataInspector::CancellationCallback& shouldCancel) {
    if (fingerprint == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal fingerprint output is null.");
        }
        return false;
    }
    const QFileInfo before(filePath);
    if (!before.exists() || !before.isFile() || before.isSymLink()
        || isFileSystemReparsePoint(before.absoluteFilePath())
        || before.size() != expectedInfo.size()
        || !sameSourcePath(before.canonicalFilePath(),
                           expectedInfo.canonicalFilePath())) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "Converted staging file changed before fingerprinting: %1")
                                .arg(filePath);
        }
        return false;
    }
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Cannot fingerprint converted staging file %1: %2")
                                .arg(filePath, file.errorString());
        }
        return false;
    }
    constexpr qint64 kHashChunkBytes = 4 * 1024 * 1024;
    QCryptographicHash hash(QCryptographicHash::Sha256);
    qint64 hashedBytes = 0;
    while (!file.atEnd()) {
        if (shouldCancel && shouldCancel()) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Converted PVD staging validation interrupted.");
            }
            return false;
        }
        const QByteArray chunk = file.read(kHashChunkBytes);
        if (chunk.isEmpty() && file.error() != QFileDevice::NoError) {
            if (errorMessage != nullptr) {
                *errorMessage = QStringLiteral("Cannot finish fingerprinting %1: %2")
                                    .arg(filePath, file.errorString());
            }
            return false;
        }
        hash.addData(chunk);
        hashedBytes += chunk.size();
        if (progress) {
            constexpr int kProgressUnits = 1000;
            const int completed = before.size() <= 0
                ? kProgressUnits
                : std::clamp(static_cast<int>(
                                 (static_cast<long double>(hashedBytes)
                                  / static_cast<long double>(before.size()))
                                 * kProgressUnits),
                             0, kProgressUnits);
            progress(QStringLiteral("Fingerprinting converted PVD staging"),
                     completed, kProgressUnits);
        }
    }
    file.close();
    const QFileInfo after(filePath);
    if (!after.exists() || !after.isFile() || after.isSymLink()
        || isFileSystemReparsePoint(after.absoluteFilePath())
        || after.size() != before.size()
        || after.lastModified() != before.lastModified()
        || !sameSourcePath(after.canonicalFilePath(), before.canonicalFilePath())) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral(
                "Converted staging file changed while fingerprinting: %1")
                                .arg(filePath);
        }
        return false;
    }
    fingerprint->path = sourcePathIdentity(filePath);
    fingerprint->size = before.size();
    fingerprint->sha256Hex = QString::fromLatin1(hash.result().toHex());
    return true;
}

FlowDataScanResult FlowDataInspector::scanConvertedPvdFolder(
    const QString& folder,
    FlowDataKind dataKind,
    bool coordinateCheck,
    const ProgressCallback& progress,
    const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataVersion = QDate::currentDate().toString(Qt::ISODate);
    result.entry.fileType = QStringLiteral("vtk series");

    auto reject = [&](const QString& message) {
        result.passed = false;
        result.errorMessage = message;
        return result;
    };
    auto canceled = [&]() {
        return shouldCancel && shouldCancel();
    };

    const QFileInfo folderInfo(folder);
    if (!folderInfo.exists() || !folderInfo.isDir() || folderInfo.isSymLink()
        || isFileSystemReparsePoint(folderInfo.absoluteFilePath())) {
        return reject(QStringLiteral(
            "Converted staging must be an existing regular, non-link directory."));
    }
    const QString stagingFolder = sourcePathIdentity(folderInfo.absoluteFilePath());
    result.sourceFolder = stagingFolder;
    result.entry.dataName = QFileInfo(stagingFolder).fileName();
    if (canceled()) {
        return reject(QStringLiteral("Converted PVD staging validation interrupted."));
    }

    QString pvdPath;
    int pvdCount = 0;
    QDirIterator pvdIterator(
        stagingFolder,
        {QStringLiteral("*.pvd")},
        QDir::AllEntries | QDir::NoDotAndDotDot | QDir::Readable
            | QDir::Hidden | QDir::System,
        QDirIterator::NoIteratorFlags);
    while (pvdIterator.hasNext()) {
        if (canceled()) {
            return reject(QStringLiteral("Converted PVD staging validation interrupted."));
        }
        pvdIterator.next();
        ++pvdCount;
        if (pvdCount > 1) {
            return reject(QStringLiteral(
                "Converted staging must contain exactly one PVD collection; found more than one."));
        }
        pvdPath = normalizedPath(pvdIterator.fileInfo().absoluteFilePath());
    }
    if (pvdCount != 1) {
        return reject(QStringLiteral(
            "Converted staging must contain exactly one PVD collection; found none."));
    }

    const QFileInfo pvdBefore(pvdPath);
    if (!pvdBefore.exists() || !pvdBefore.isFile() || pvdBefore.isSymLink()
        || isFileSystemReparsePoint(pvdBefore.absoluteFilePath())
        || !pvdBefore.isReadable()
        || !sameSourcePath(pvdBefore.absolutePath(), stagingFolder)
        || pvdBefore.size() <= 0
        || pvdBefore.size()
            > Streamcenter::FlowDataProtocolLimits::MaximumCollectionMetadataBytes) {
        return reject(QStringLiteral(
            "Converted PVD collection is missing, linked, unreadable, escaped, empty, or exceeds the bounded metadata size."));
    }
    QFile pvdFile(pvdPath);
    if (!pvdFile.open(QIODevice::ReadOnly)) {
        return reject(QStringLiteral("Cannot read converted PVD collection: %1")
                          .arg(pvdFile.errorString()));
    }
    const QByteArray pvdPayload = pvdFile.read(
        Streamcenter::FlowDataProtocolLimits::MaximumCollectionMetadataBytes + 1);
    if (pvdPayload.size()
            > Streamcenter::FlowDataProtocolLimits::MaximumCollectionMetadataBytes
        || !pvdFile.atEnd()) {
        return reject(QStringLiteral(
            "Converted PVD collection exceeds the bounded metadata size."));
    }
    pvdFile.close();
    const QFileInfo pvdAfterRead(pvdPath);
    if (!pvdAfterRead.exists() || !pvdAfterRead.isFile() || pvdAfterRead.isSymLink()
        || isFileSystemReparsePoint(pvdAfterRead.absoluteFilePath())
        || pvdAfterRead.size() != pvdBefore.size()
        || pvdAfterRead.lastModified() != pvdBefore.lastModified()
        || !sameSourcePath(pvdAfterRead.canonicalFilePath(),
                           pvdBefore.canonicalFilePath())) {
        return reject(QStringLiteral(
            "Converted PVD collection changed while it was being read."));
    }

    QVector<double> timeValues;
    QStringList referencedFiles;
    QSet<QString> referencedIdentities;
    QStringList elementStack;
    bool rootSeen = false;
    bool rootClosed = false;
    bool collectionSeen = false;
    QXmlStreamReader xml(pvdPayload);
    while (!xml.atEnd()) {
        if (canceled()) {
            return reject(QStringLiteral("Converted PVD staging validation interrupted."));
        }
        xml.readNext();
        if (xml.isStartElement()) {
            const QString elementName = xml.name().toString();
            if (!xml.namespaceUri().isEmpty()) {
                return reject(QStringLiteral(
                    "Converted PVD must not use namespaced structural elements."));
            }
            if (elementStack.isEmpty()) {
                if (rootSeen || rootClosed || elementName != QStringLiteral("VTKFile")
                    || xml.attributes().value(QStringLiteral("type"))
                        != QStringLiteral("Collection")) {
                    return reject(QStringLiteral(
                        "Converted PVD root must be exactly one VTKFile of type Collection."));
                }
                rootSeen = true;
                elementStack.append(elementName);
                continue;
            }
            if (elementStack.size() == 1) {
                if (elementStack.constLast() != QStringLiteral("VTKFile")
                    || collectionSeen
                    || elementName != QStringLiteral("Collection")) {
                    return reject(QStringLiteral(
                        "Converted PVD VTKFile must contain exactly one direct Collection element."));
                }
                collectionSeen = true;
                elementStack.append(elementName);
                continue;
            }
            if (elementStack.size() != 2
                || elementStack.at(0) != QStringLiteral("VTKFile")
                || elementStack.at(1) != QStringLiteral("Collection")
                || elementName != QStringLiteral("DataSet")) {
                return reject(QStringLiteral(
                    "Converted PVD DataSet elements must be direct children of its sole Collection."));
            }
            if (timeValues.size()
                >= Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount) {
                return reject(QStringLiteral("Converted PVD exceeds the %1-frame limit.")
                                  .arg(Streamcenter::FlowDataProtocolLimits::MaximumResultFileCount));
            }
            const QString timeText =
                xml.attributes().value(QStringLiteral("timestep")).toString();
            const QString relativePath =
                xml.attributes().value(QStringLiteral("file")).toString();
            bool timeOk = false;
            const double timeValue = timeText.toDouble(&timeOk);
            const QString normalizedRelativePath =
                QDir::fromNativeSeparators(relativePath);
            const QString cleanRelativePath = QDir::cleanPath(normalizedRelativePath);
            if (!timeOk || !std::isfinite(timeValue) || relativePath.isEmpty()
                || relativePath.toUtf8().size()
                    > Streamcenter::FlowDataProtocolLimits::MaximumStringBytes
                || QDir::isAbsolutePath(relativePath)
                || cleanRelativePath != normalizedRelativePath
                || QFileInfo(cleanRelativePath).path() != QStringLiteral(".")) {
                return reject(QStringLiteral(
                    "Converted PVD contains an invalid timestep or non-direct file reference."));
            }
            if (QFileInfo(cleanRelativePath).suffix().compare(
                    QStringLiteral("vtu"), Qt::CaseInsensitive) != 0) {
                return reject(QStringLiteral(
                    "Converted PVD may reference only converter-produced VTU frames: %1")
                                  .arg(relativePath));
            }
            const QFileInfo frameInfo(
                QDir(stagingFolder).absoluteFilePath(cleanRelativePath));
            if (!frameInfo.exists() || !frameInfo.isFile() || frameInfo.isSymLink()
                || isFileSystemReparsePoint(frameInfo.absoluteFilePath())
                || !frameInfo.isReadable()
                || !sameSourcePath(frameInfo.absolutePath(), stagingFolder)) {
                return reject(QStringLiteral(
                    "Converted PVD references a missing, linked, unreadable, or escaped frame: %1")
                                  .arg(relativePath));
            }
            const QString absoluteFramePath =
                sourcePathIdentity(frameInfo.absoluteFilePath());
            const QString frameIdentity = sourcePathIdentityKey(absoluteFramePath);
            if (referencedIdentities.contains(frameIdentity)) {
                return reject(QStringLiteral(
                    "Converted PVD references a frame more than once: %1")
                                  .arg(relativePath));
            }
            referencedIdentities.insert(frameIdentity);
            timeValues.append(timeValue);
            referencedFiles.append(absoluteFramePath);
            elementStack.append(elementName);
            continue;
        }
        if (xml.isEndElement()) {
            const QString elementName = xml.name().toString();
            if (elementStack.isEmpty() || elementStack.constLast() != elementName) {
                return reject(QStringLiteral(
                    "Converted PVD XML element nesting is inconsistent."));
            }
            elementStack.removeLast();
            if (elementStack.isEmpty()) {
                rootClosed = true;
            }
            continue;
        }
        if (xml.isCharacters() && !xml.isWhitespace()) {
            return reject(QStringLiteral(
                "Converted PVD structural elements may contain only whitespace text."));
        }
    }
    if (xml.hasError()) {
        return reject(QStringLiteral("Converted PVD XML is malformed: %1")
                          .arg(xml.errorString()));
    }
    if (!rootSeen || !rootClosed || !collectionSeen || !elementStack.isEmpty()
        || timeValues.isEmpty()) {
        return reject(QStringLiteral(
            "Converted PVD is not a complete, non-empty VTK Collection."));
    }

    const QString pvdIdentity = sourcePathIdentityKey(pvdPath);
    QSet<QString> expectedStagingIdentities = referencedIdentities;
    expectedStagingIdentities.insert(pvdIdentity);
    auto validateStagingEntries = [&]() {
        QSet<QString> seen;
        QDirIterator iterator(
            stagingFolder,
            QDir::AllEntries | QDir::NoDotAndDotDot | QDir::Hidden | QDir::System,
            QDirIterator::NoIteratorFlags);
        while (iterator.hasNext()) {
            if (canceled()) {
                result.errorMessage =
                    QStringLiteral("Converted PVD staging validation interrupted.");
                return false;
            }
            iterator.next();
            const QFileInfo entry = iterator.fileInfo();
            if (seen.size() >= expectedStagingIdentities.size()) {
                result.errorMessage = QStringLiteral(
                    "Converted staging contains more entries than its bounded PVD manifest permits.");
                return false;
            }
            const QString identity = sourcePathIdentityKey(entry.absoluteFilePath());
            if (!entry.isFile() || entry.isSymLink()
                || isFileSystemReparsePoint(entry.absoluteFilePath())
                || !expectedStagingIdentities.contains(identity)
                || seen.contains(identity)) {
                result.errorMessage = QStringLiteral(
                    "Converted staging contains an extra, duplicate, linked, or non-file entry outside the PVD manifest: %1")
                                          .arg(entry.absoluteFilePath());
                return false;
            }
            seen.insert(identity);
        }
        if (seen != expectedStagingIdentities) {
            result.errorMessage = QStringLiteral(
                "Converted staging does not contain exactly the files declared by its PVD manifest.");
            return false;
        }
        return true;
    };
    if (!validateStagingEntries()) {
        return reject(result.errorMessage);
    }

    double dt = 1.0;
    if (timeValues.size() > 1) {
        dt = timeValues.at(1) - timeValues.at(0);
        if (!std::isfinite(dt) || dt <= 0.0) {
            return reject(QStringLiteral(
                "Converted PVD timesteps must be strictly increasing."));
        }
        for (qsizetype index = 2; index < timeValues.size(); ++index) {
            const double current = timeValues.at(index);
            const double previous = timeValues.at(index - 1);
            const double expected = std::fma(static_cast<double>(index),
                                             dt, timeValues.first());
            const double magnitude = std::max(std::abs(current), std::abs(expected));
            const double tolerance = std::abs(dt)
                    * static_cast<double>(index) * 1.0e-10
                + std::numeric_limits<double>::epsilon() * magnitude * 8.0
                + std::numeric_limits<double>::denorm_min() * 16.0;
            if (!std::isfinite(current) || !(current > previous)
                || !std::isfinite(expected)
                || std::abs(current - expected) > tolerance) {
                return reject(QStringLiteral(
                    "Converted PVD uses a non-uniform time axis that cannot be represented by one dt value."));
            }
        }
    }

    auto effectiveDimension = [](const VtkDataSummary& summary) {
        return summary.spatialDimension;
    };
    auto sameCoordinates = [](const VtkDataSummary& lhs,
                              const VtkDataSummary& rhs) {
        return coordinatesEquivalent(lhs, rhs);
    };

    VtkDataSummary firstFrame;
    int firstDimension = 0;
    QVector<QFileInfo> validatedFrameInfos;
    validatedFrameInfos.reserve(referencedFiles.size());
    QVector<FlowDataSourceFingerprint> preParseFrameFingerprints;
    preParseFrameFingerprints.reserve(referencedFiles.size());
    for (qsizetype index = 0; index < referencedFiles.size(); ++index) {
        if (canceled()) {
            return reject(QStringLiteral("Converted PVD staging validation interrupted."));
        }
        if (progress) {
            progress(QStringLiteral("Validating converted PVD frames"),
                     static_cast<int>(index),
                     static_cast<int>(referencedFiles.size()));
        }
        const QFileInfo frameBefore(referencedFiles.at(index));
        FlowDataSourceFingerprint preParseFingerprint;
        QString fingerprintError;
        if (!fingerprintStableRegularFile(referencedFiles.at(index), frameBefore,
                                          &preParseFingerprint,
                                          &fingerprintError,
                                          progress, shouldCancel)) {
            return reject(fingerprintError);
        }
        QString frameError;
        VtkDataSummary frame;
        if (!summarizeAnyVtk(referencedFiles.at(index), &frame, &frameError)) {
            return reject(frameError);
        }
        const QFileInfo frameAfter(referencedFiles.at(index));
        if (!frameAfter.exists() || !frameAfter.isFile() || frameAfter.isSymLink()
            || isFileSystemReparsePoint(frameAfter.absoluteFilePath())
            || frameAfter.size() != frameBefore.size()
            || frameAfter.lastModified() != frameBefore.lastModified()
            || !sameSourcePath(frameAfter.canonicalFilePath(),
                               frameBefore.canonicalFilePath())) {
            return reject(QStringLiteral(
                "Converted VTU frame changed while it was being validated: %1")
                              .arg(referencedFiles.at(index)));
        }
        const int frameDimension = effectiveDimension(frame);
        if (dataKind == FlowDataKind::Resampled
            && (!frame.uniformStructured || !frame.hasDims)) {
            return reject(QStringLiteral(
                "Every converted resampled frame must be a uniform structured VTK dataset: %1")
                              .arg(referencedFiles.at(index)));
        }
        if (index == 0) {
            firstFrame = frame;
            firstDimension = frameDimension;
        } else {
            if (!sameCandidateFieldSet(firstFrame.termCandidates,
                                       frame.termCandidates)
                || !sameArraySchema(firstFrame.arrays, frame.arrays,
                                    dataKind == FlowDataKind::Resampled)) {
                return reject(QStringLiteral(
                    "Converted PVD field/component schema differs from the first frame: %1")
                                  .arg(referencedFiles.at(index)));
            }
            if (frameDimension != firstDimension) {
                return reject(QStringLiteral(
                    "Converted PVD spatial dimension differs from the first frame: %1")
                                  .arg(referencedFiles.at(index)));
            }
            if (dataKind == FlowDataKind::Resampled
                && (!firstFrame.hasDims || !frame.hasDims
                    || firstFrame.dims[0] != frame.dims[0]
                    || firstFrame.dims[1] != frame.dims[1]
                    || firstFrame.dims[2] != frame.dims[2])) {
                return reject(QStringLiteral(
                    "Converted resampled PVD grid dimensions differ from the first frame: %1")
                                  .arg(referencedFiles.at(index)));
            }
            if (coordinateCheck && !sameCoordinates(firstFrame, frame)
                && result.warningMessage.isEmpty()) {
                result.warningMessage = QStringLiteral(
                    "Converted PVD contains changed mesh coordinates; dimensions and all frame field schemas were still validated.");
            }
        }
        validatedFrameInfos.append(frameAfter);
        preParseFrameFingerprints.append(preParseFingerprint);
        if (progress) {
            progress(QStringLiteral("Validating converted PVD frames"),
                     static_cast<int>(index + 1),
                     static_cast<int>(referencedFiles.size()));
        }
        if (canceled()) {
            return reject(QStringLiteral("Converted PVD staging validation interrupted."));
        }
    }
    if (firstFrame.termCandidates.isEmpty()) {
        return reject(QStringLiteral(
            "Converted PVD frames do not contain numeric point or cell arrays."));
    }

    auto fingerprintFile = [&](const QString& filePath,
                               const QFileInfo& expectedInfo,
                               FlowDataSourceFingerprint* fingerprint) {
        if (fingerprint == nullptr) {
            result.errorMessage = QStringLiteral(
                "Internal error: converted staging fingerprint output is null.");
            return false;
        }
        const QFileInfo before(filePath);
        if (!before.exists() || !before.isFile() || before.isSymLink()
            || isFileSystemReparsePoint(before.absoluteFilePath())
            || before.size() != expectedInfo.size()
            || before.lastModified() != expectedInfo.lastModified()
            || !sameSourcePath(before.canonicalFilePath(),
                               expectedInfo.canonicalFilePath())) {
            result.errorMessage = QStringLiteral(
                "Converted staging file changed before it could be fingerprinted: %1")
                                      .arg(filePath);
            return false;
        }
        QFile file(filePath);
        if (!file.open(QIODevice::ReadOnly)) {
            result.errorMessage = QStringLiteral(
                "Cannot fingerprint converted staging file %1: %2")
                                      .arg(filePath, file.errorString());
            return false;
        }
        QCryptographicHash hash(QCryptographicHash::Sha256);
        constexpr qint64 kHashChunkBytes = 4 * 1024 * 1024;
        qint64 hashedBytes = 0;
        while (!file.atEnd()) {
            if (canceled()) {
                result.errorMessage =
                    QStringLiteral("Converted PVD staging validation interrupted.");
                return false;
            }
            const QByteArray chunk = file.read(kHashChunkBytes);
            if (chunk.isEmpty() && file.error() != QFileDevice::NoError) {
                result.errorMessage = QStringLiteral(
                    "Cannot finish fingerprinting converted staging file %1: %2")
                                          .arg(filePath, file.errorString());
                return false;
            }
            hash.addData(chunk);
            hashedBytes += chunk.size();
            if (progress) {
                constexpr int kProgressUnits = 1000;
                const int completedUnits = before.size() <= 0
                    ? kProgressUnits
                    : std::clamp(
                          static_cast<int>(
                              (static_cast<long double>(hashedBytes)
                               / static_cast<long double>(before.size()))
                              * kProgressUnits),
                          0,
                          kProgressUnits);
                // Besides user-visible progress, every bounded read emits a
                // protocol heartbeat.  This prevents a healthy large-file
                // SHA-256 pass from tripping the parent helper's inactivity
                // timeout while keeping all heavy I/O outside the GUI process.
                progress(QStringLiteral("Fingerprinting converted PVD staging"),
                         completedUnits,
                         kProgressUnits);
            }
        }
        file.close();
        const QFileInfo after(filePath);
        if (!after.exists() || !after.isFile() || after.isSymLink()
            || isFileSystemReparsePoint(after.absoluteFilePath())
            || after.size() != before.size()
            || after.lastModified() != before.lastModified()
            || !sameSourcePath(after.canonicalFilePath(), before.canonicalFilePath())) {
            result.errorMessage = QStringLiteral(
                "Converted staging file changed while it was being fingerprinted: %1")
                                      .arg(filePath);
            return false;
        }
        fingerprint->path = sourcePathIdentity(filePath);
        fingerprint->size = before.size();
        fingerprint->sha256Hex = QString::fromLatin1(hash.result().toHex());
        return true;
    };

    result.sourceInventory.clear();
    result.sourceInventory.reserve(referencedFiles.size() + 1);
    FlowDataSourceFingerprint pvdFingerprint;
    if (!fingerprintFile(pvdPath, pvdAfterRead, &pvdFingerprint)) {
        return reject(result.errorMessage);
    }
    const QString parsedPvdHash = QString::fromLatin1(
        QCryptographicHash::hash(pvdPayload, QCryptographicHash::Sha256).toHex());
    if (pvdFingerprint.sha256Hex != parsedPvdHash) {
        return reject(QStringLiteral(
            "Converted PVD collection changed after it was parsed."));
    }
    result.sourceInventory.append(pvdFingerprint);
    for (qsizetype index = 0; index < referencedFiles.size(); ++index) {
        FlowDataSourceFingerprint fingerprint;
        if (!fingerprintFile(referencedFiles.at(index),
                             validatedFrameInfos.at(index),
                             &fingerprint)) {
            return reject(result.errorMessage);
        }
        const FlowDataSourceFingerprint& beforeParse =
            preParseFrameFingerprints.at(index);
        if (fingerprint.path != beforeParse.path
            || fingerprint.size != beforeParse.size
            || fingerprint.sha256Hex != beforeParse.sha256Hex) {
            return reject(QStringLiteral(
                "Converted VTU frame content changed across native VTK parsing: %1")
                              .arg(referencedFiles.at(index)));
        }
        result.sourceInventory.append(fingerprint);
        if (progress) {
            progress(QStringLiteral("Fingerprinting converted PVD staging"),
                     static_cast<int>(index + 1),
                     static_cast<int>(referencedFiles.size()));
        }
    }
    std::sort(result.sourceInventory.begin(), result.sourceInventory.end(),
              [](const FlowDataSourceFingerprint& lhs,
                 const FlowDataSourceFingerprint& rhs) {
                  return sourcePathIdentityKey(lhs.path)
                      < sourcePathIdentityKey(rhs.path);
              });
    if (!validateStagingEntries()) {
        return reject(result.errorMessage);
    }

    result.files = referencedFiles;
    result.termCandidates = firstFrame.termCandidates;
    result.entry.termCandidates = result.termCandidates;
    inferTerms(firstFrame.arrays, firstDimension, &result.entry);
    result.entry.xRange = formatRange(firstFrame.bounds[0], firstFrame.bounds[1]);
    result.entry.yRange = formatRange(firstFrame.bounds[2], firstFrame.bounds[3]);
    result.entry.zRange = formatRange(firstFrame.bounds[4], firstFrame.bounds[5]);
    result.entry.dimension = QString::number(firstDimension);
    result.entry.meshCount = dataKind == FlowDataKind::Resampled
        ? QStringLiteral("%1 %2 %3")
              .arg(std::max(1, firstFrame.dims[0]))
              .arg(std::max(1, firstFrame.dims[1]))
              .arg(std::max(1, firstFrame.dims[2]))
        : QString::number(referencedFiles.size());
    result.entry.dt = formatDouble(dt);
    result.entry.checkStatus = QStringLiteral("passed");
    result.passed = true;
    return result;
}

FlowDataScanResult FlowDataInspector::scanNek5000File(const QString& filePath,
                                                      FlowDataKind dataKind,
                                                      bool coordinateCheck,
                                                      const ProgressCallback& progress,
                                                      const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    const QString normalizedFilePath = normalizedPath(filePath);
    result.sourceFolder = normalizedFilePath;
    result.files = {normalizedFilePath};
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataName = QFileInfo(filePath).completeBaseName();
    result.entry.dataVersion = QDate::currentDate().toString(Qt::ISODate);
    result.entry.fileType = QStringLiteral("nek5000 series");

    auto canceled = [&]() {
        if (shouldCancel && shouldCancel()) {
            result.errorMessage = QStringLiteral("Nek5000 validation interrupted.");
            return true;
        }
        return false;
    };
    if (canceled()) {
        return result;
    }

    const QFileInfo fileInfo(filePath);
    if (!fileInfo.exists() || !fileInfo.isFile()) {
        result.errorMessage = QString("Nek5000 descriptor does not exist: %1").arg(filePath);
        return result;
    }
    const QString suffix = fileInfo.suffix();
    if (suffix.compare(QStringLiteral("nek5000"), Qt::CaseInsensitive) != 0 &&
        suffix.compare(QStringLiteral("nek"), Qt::CaseInsensitive) != 0) {
        result.errorMessage = QString("Nek5000 import expects a *.nek5000 or *.nek descriptor file: %1")
                                  .arg(filePath);
        return result;
    }

    QString readError;
    std::vector<streamcenterplus::io::DataSetFrame> frames;
    if (!enumerateNek5000Frames(normalizedFilePath, &frames, &readError)) {
        result.errorMessage = readError;
        return result;
    }
    QVector<double> timeValues;
    timeValues.reserve(static_cast<qsizetype>(frames.size()));
    // Do not sort or deduplicate: source index is part of a Native Nek5000
    // frame's identity, even when two snapshots advertise the same time.
    for (const streamcenterplus::io::DataSetFrame& frame : frames) {
        timeValues.append(frame.time);
    }

    const int progressTotal = coordinateCheck
        ? static_cast<int>(std::min<std::size_t>(
              frames.size() + 1, static_cast<std::size_t>(std::numeric_limits<int>::max())))
        : 2;
    if (progress) {
        progress(QStringLiteral("Validating Nek5000 flow data"), 1, progressTotal);
    }
    if (canceled()) {
        return result;
    }

    if (timeValues.size() >= 2) {
        const double dt = timeValues.at(1) - timeValues.at(0);
        bool uniform = std::isfinite(dt) && dt > 0.0;
        for (int index = 2; index < timeValues.size(); ++index) {
            if (std::abs((timeValues.at(index) - timeValues.at(index - 1)) - dt) > 1.0e-10) {
                uniform = false;
                break;
            }
        }
        if (uniform) {
            result.entry.dt = formatDouble(dt);
        }
    }

    vtkSmartPointer<vtkDataSet> first = readCanonicalNek5000Frame(frames.front(), &readError);
    if (first == nullptr) {
        result.errorMessage = readError;
        return result;
    }

    int dims[3] = {0, 0, 0};
    if (dataKind == FlowDataKind::Resampled && !isUniformStructured(first, dims)) {
        result.errorMessage = "Resampled flow data must be a uniform structured VTK dataset.";
        return result;
    }
    if (dataKind != FlowDataKind::Resampled) {
        dimensionsForDataSet(first, dims);
    }

    double bounds[6] = {};
    if (!canonicalNek5000Bounds(first, bounds, &readError)) {
        result.errorMessage = readError;
        return result;
    }
    double nodalBounds[6] = {};
    first->GetBounds(nodalBounds);
    if (!coordinateBoundsEquivalent(bounds, nodalBounds)) {
        result.warningMessage = QStringLiteral(
            "Native Nek5000 coordinate ranges conservatively enclose the complete "
            "high-order spectral mapping; they are not limited to the raw GLL-node bounds.");
    }
    result.entry.xRange = formatRange(bounds[0], bounds[1]);
    result.entry.yRange = formatRange(bounds[2], bounds[3]);
    result.entry.zRange = formatRange(bounds[4], bounds[5]);
    const int dimension = spatialDimensionForDataSet(first, dims);
    result.entry.dimension = QString::number(dimension);

    if (dataKind == FlowDataKind::Resampled) {
        result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                     .arg(std::max(1, dims[0]))
                                     .arg(std::max(1, dims[1]))
                                     .arg(std::max(1, dims[2]));
    } else {
        result.entry.meshCount = QString::number(static_cast<qlonglong>(frames.size()));
    }

    QList<ArrayInfo> arrays;
    collectArrays(first->GetPointData(), ArrayInfo::Association::Point,
                  &arrays, &result.termCandidates);
    collectArrays(first->GetCellData(), ArrayInfo::Association::Cell,
                  &arrays, &result.termCandidates);
    result.termCandidates.sort(Qt::CaseInsensitive);
    inferTerms(arrays, dimension, &result.entry);
    result.entry.termCandidates = result.termCandidates;

    if (result.termCandidates.isEmpty()) {
        result.errorMessage = "The Nek5000 snapshot does not contain numeric point or cell arrays.";
        return result;
    }

    if (coordinateCheck) {
        for (std::size_t index = 1; index < frames.size(); ++index) {
            if (canceled()) {
                return result;
            }
            vtkSmartPointer<vtkDataSet> frame = readCanonicalNek5000Frame(
                frames[index], &readError);
            if (frame == nullptr) {
                result.errorMessage = readError;
                return result;
            }
            if (!streamcenterplus::io::HasSameNek5000SpectralGeometry(first, frame)) {
                result.errorMessage = QString(
                    "Nek5000 canonical element-local GLL geometry differs from the first frame: %1")
                                          .arg(qStringFromFsPath(frames[index].path));
                return result;
            }
            if (progress) {
                progress(QStringLiteral("Validating Nek5000 flow data"),
                         static_cast<int>(index + 2), progressTotal);
            }
            if (canceled()) {
                return result;
            }
        }
    }

    if (progress) {
        progress(QStringLiteral("Validating Nek5000 flow data"), progressTotal, progressTotal);
    }
    if (canceled()) {
        return result;
    }

    result.entry.checkStatus = QStringLiteral("passed");
    result.passed = true;
    return result;
}

FlowDataScanResult FlowDataInspector::scanOpenFoamCase(const QString& casePath,
                                                       const QString& caseType,
                                                       FlowDataKind dataKind,
                                                       bool coordinateCheck,
                                                       bool startEndFieldsOnly,
                                                       const ProgressCallback& progress,
                                                       const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    result.sourceFolder = normalizedPath(casePath);
    result.files = {normalizedPath(casePath)};
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataName = QFileInfo(casePath).completeBaseName();
    result.entry.dataVersion = QDate::currentDate().toString(Qt::ISODate);
    result.entry.fileType = openFoamFileType(caseType);

    auto canceled = [&]() {
        if (shouldCancel && shouldCancel()) {
            result.errorMessage = QStringLiteral("OpenFOAM validation interrupted.");
            return true;
        }
        return false;
    };
    if (canceled()) {
        return result;
    }

    double writeInterval = 0.0;
    if (!readOpenFoamWriteInterval(casePath, &writeInterval, &result.errorMessage)) {
        return result;
    }
    result.entry.dt = formatDouble(writeInterval);

    vtkSmartPointer<vtkOpenFOAMReader> reader = createOpenFoamReader(caseType);
    if (!prepareOpenFoamReader(reader, casePath, caseType, &result.errorMessage)) {
        return result;
    }
    if (canceled()) {
        return result;
    }

    const QVector<double> timeValues = openFoamTimeValues(reader);
    if (timeValues.isEmpty()) {
        result.errorMessage = QStringLiteral("The OpenFOAM case does not contain any time steps.");
        return result;
    }

    const bool validateEveryFrame = coordinateCheck || startEndFieldsOnly;
    if (progress) {
        progress(QStringLiteral("Validating OpenFOAM case"), 1, validateEveryFrame ? timeValues.size() : 2);
    }
    if (canceled()) {
        return result;
    }

    QString readError;
    vtkSmartPointer<vtkDataSet> first = readOpenFoamVolumeAtTime(reader, timeValues.first(), &readError);
    if (first == nullptr) {
        result.errorMessage = readError;
        return result;
    }

    int dims[3] = {0, 0, 0};
    if (dataKind == FlowDataKind::Resampled && !isUniformStructured(first, dims)) {
        result.errorMessage = "Resampled flow data must be a uniform structured VTK dataset.";
        return result;
    }
    if (dataKind != FlowDataKind::Resampled) {
        dimensionsForDataSet(first, dims);
    }

    double bounds[6] = {};
    first->GetBounds(bounds);
    result.entry.xRange = formatRange(bounds[0], bounds[1]);
    result.entry.yRange = formatRange(bounds[2], bounds[3]);
    result.entry.zRange = formatRange(bounds[4], bounds[5]);
    const int dimension = spatialDimensionForDataSet(first, dims);
    result.entry.dimension = QString::number(dimension);

    if (dataKind == FlowDataKind::Resampled) {
        result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                     .arg(std::max(1, dims[0]))
                                     .arg(std::max(1, dims[1]))
                                     .arg(std::max(1, dims[2]));
    } else {
        result.entry.meshCount = QString::number(timeValues.size());
    }

    QList<ArrayInfo> arrays;
    collectArrays(first->GetPointData(), ArrayInfo::Association::Point,
                  &arrays, &result.termCandidates);
    collectArrays(first->GetCellData(), ArrayInfo::Association::Cell,
                  &arrays, &result.termCandidates);
    QStringList commonFields = actualArrayNames(first);
    if (commonFields.isEmpty()) {
        result.errorMessage = "The first OpenFOAM frame does not contain numeric point or cell arrays.";
        return result;
    }

    if (validateEveryFrame) {
        for (int index = 1; index < timeValues.size(); ++index) {
            if (canceled()) {
                return result;
            }
            QString frameError;
            vtkSmartPointer<vtkDataSet> frame = readOpenFoamVolumeAtTime(reader, timeValues.at(index), &frameError);
            if (frame == nullptr) {
                result.errorMessage = frameError;
                return result;
            }
            if (coordinateCheck && !coordinatesEquivalent(first, frame)) {
                result.errorMessage = QString("OpenFOAM frame grid coordinates differ from the first frame at time %1.")
                                          .arg(formatDouble(timeValues.at(index)));
                return result;
            }
            if (startEndFieldsOnly) {
                commonFields = intersectFieldNames(commonFields, actualArrayNames(frame));
                if (commonFields.isEmpty()) {
                    result.errorMessage = QStringLiteral("No numeric OpenFOAM fields are present in every time step.");
                    return result;
                }
            }
            if (progress) {
                progress(QStringLiteral("Validating OpenFOAM case"), index + 1, timeValues.size());
            }
            if (canceled()) {
                return result;
            }
        }
    } else if (progress) {
        progress(QStringLiteral("Validating OpenFOAM case"), 2, 2);
        if (canceled()) {
            return result;
        }
    }

    if (startEndFieldsOnly) {
        arrays = filterArrayInfos(arrays, commonFields);
    }
    result.termCandidates = candidateListForArrays(arrays);
    inferTerms(arrays, dimension, &result.entry);
    result.entry.termCandidates = result.termCandidates;

    if (result.termCandidates.isEmpty()) {
        result.errorMessage = startEndFieldsOnly
            ? QStringLiteral("No numeric OpenFOAM fields are present in every time step.")
            : QStringLiteral("The first OpenFOAM frame does not contain numeric point or cell arrays.");
        return result;
    }

    result.entry.checkStatus = QStringLiteral("passed");
    result.passed = true;
    return result;
}

FlowDataScanResult FlowDataInspector::scanEnSightGoldCase(const QString& casePath,
                                                          FlowDataKind dataKind,
                                                          bool coordinateCheck,
                                                          bool commonFieldsOnly,
                                                          const ProgressCallback& progress,
                                                          const CancellationCallback& shouldCancel) {
    return scanEnSightGoldCases({casePath}, dataKind, coordinateCheck, commonFieldsOnly, progress, shouldCancel);
}

FlowDataScanResult FlowDataInspector::scanEnSightGoldCases(const QStringList& casePaths,
                                                           FlowDataKind dataKind,
                                                           bool coordinateCheck,
                                                           bool commonFieldsOnly,
                                                           const ProgressCallback& progress,
                                                           const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    result.sourceFolder = casePaths.size() == 1 ? normalizedPath(casePaths.first()) : QString();
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataName = casePaths.size() == 1 ? QFileInfo(casePaths.first()).completeBaseName() : QStringLiteral("ensight_gold");
    result.entry.dataVersion = QDate::currentDate().toString(Qt::ISODate);
    result.entry.fileType = enSightGoldFileType();

    auto canceled = [&]() {
        if (shouldCancel && shouldCancel()) {
            result.errorMessage = QStringLiteral("EnSight Gold validation interrupted.");
            return true;
        }
        return false;
    };
    if (canceled()) {
        return result;
    }

    QVector<EnSightGoldCaseSource> cases;
    QVector<EnSightGoldFrameSource> frames;
    if (!loadEnSightGoldCaseSources(casePaths, &cases, &frames, &result.errorMessage)) {
        return result;
    }
    if (!buildEnSightGoldInventoryFromCases(cases,
                                            &result.sourceInventory,
                                            &result.errorMessage,
                                            progress,
                                            shouldCancel)) {
        return result;
    }
    result.files.reserve(cases.size());
    for (const EnSightGoldCaseSource& source : std::as_const(cases)) {
        result.files.append(source.path);
    }
    if (cases.size() > 1) {
        const QString commonParent = normalizedPath(QFileInfo(cases.first().path).absolutePath());
        const bool allShareParent = std::all_of(
            cases.cbegin(), cases.cend(), [&](const EnSightGoldCaseSource& source) {
                return sameSourcePath(QFileInfo(source.path).absolutePath(), commonParent);
            });
        result.sourceFolder = allShareParent ? commonParent : QString();
    }
    if (canceled()) {
        return result;
    }

    if (frames.isEmpty()) {
        result.errorMessage = QStringLiteral("The EnSight Gold case does not contain any readable frames.");
        return result;
    }

    QVector<double> timeValues;
    timeValues.reserve(frames.size());
    for (const EnSightGoldFrameSource& frame : std::as_const(frames)) {
        timeValues.append(frame.timeValue);
    }
    double dt = 0.0;
    if (!uniformTimeInterval(timeValues, &dt, &result.errorMessage)) {
        return result;
    }
    result.entry.dt = formatDouble(dt);

    const bool validateEveryFrame = coordinateCheck || commonFieldsOnly;
    if (progress) {
        progress(QStringLiteral("Validating EnSight Gold case"), 1, validateEveryFrame ? frames.size() : 2);
    }
    if (canceled()) {
        return result;
    }

    if (anyEnSightGoldSourceLooksLarge(cases)) {
        QList<ArrayInfo> arrays = cases.first().arrays;
        if (commonFieldsOnly) {
            QStringList commonFields = arrayInfoNames(arrays);
            for (int index = 1; index < cases.size(); ++index) {
                commonFields = intersectFieldNames(commonFields, arrayInfoNames(cases.at(index).arrays));
                if (commonFields.isEmpty()) {
                    result.errorMessage = QStringLiteral("No numeric EnSight Gold variables are declared in every selected *.case file.");
                    return result;
                }
            }
            arrays = filterArrayInfos(arrays, commonFields);
        }

        result.termCandidates = candidateListForArrays(arrays);
        result.entry.termCandidates = result.termCandidates;
        result.entry.meshCount = QString::number(frames.size());
        // Native geometry was intentionally not opened, so dimension, ranges,
        // and inferred vector components are unknown.  Leaving them empty
        // prevents a metadata-only 2D source from being falsely pinned to 3D;
        // the isolated converted-PVD scan supplies authoritative values.
        result.entry.dimension.clear();
        result.entry.xRange.clear();
        result.entry.yRange.clear();
        result.entry.zRange.clear();
        if (result.termCandidates.isEmpty()) {
            result.errorMessage = commonFieldsOnly
                ? QStringLiteral("No numeric EnSight Gold variables are declared in every selected *.case file.")
                : QStringLiteral("The EnSight Gold case does not declare numeric variables.");
            return result;
        }
        result.warningMessage = coordinateCheck
            ? QStringLiteral("EnSight Gold validation used *.case metadata only because the referenced native files are large. Coordinate ranges were not loaded during validation.")
            : QStringLiteral("EnSight Gold validation used *.case metadata only because the referenced native files are large.");
        if (progress) {
            progress(QStringLiteral("Validating EnSight Gold case"), 2, 2);
            if (canceled()) {
                return result;
            }
        }
        result.entry.checkStatus = QStringLiteral("passed");
        result.passed = true;
        return result;
    }

    QList<ArrayInfo> arrays;
    QStringList commonFields;
    double firstBounds[6] = {};
    QString firstBoundsText;
    QByteArray firstCoordinateSignature;
    int dimension = 3;
    {
        const EnSightGoldFrameSource& firstFrame = frames.first();
        const QString firstCasePath = cases.at(firstFrame.caseIndex).path;
        QString readError;
        vtkSmartPointer<vtkDataSet> first = readEnSightGoldVolumeFromCaseAtTime(firstCasePath,
                                                                                firstFrame.timeValue,
                                                                                &readError);
        if (first == nullptr) {
            result.errorMessage = readError;
            return result;
        }

        int dims[3] = {0, 0, 0};
        if (dataKind == FlowDataKind::Resampled && !isUniformStructured(first, dims)) {
            result.errorMessage = "Resampled flow data must be a uniform structured VTK dataset.";
            return result;
        }
        if (dataKind != FlowDataKind::Resampled) {
            dimensionsForDataSet(first, dims);
        }

        first->GetBounds(firstBounds);
        firstCoordinateSignature = coordinateSignature(first);
        firstBoundsText = boundsSummary(firstBounds);
        result.entry.xRange = formatRange(firstBounds[0], firstBounds[1]);
        result.entry.yRange = formatRange(firstBounds[2], firstBounds[3]);
        result.entry.zRange = formatRange(firstBounds[4], firstBounds[5]);
        dimension = spatialDimensionForDataSet(first, dims);
        result.entry.dimension = QString::number(dimension);

        if (dataKind == FlowDataKind::Resampled) {
            result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                         .arg(std::max(1, dims[0]))
                                         .arg(std::max(1, dims[1]))
                                         .arg(std::max(1, dims[2]));
        } else {
            result.entry.meshCount = QString::number(frames.size());
        }

        collectArrays(first->GetPointData(), ArrayInfo::Association::Point,
                      &arrays, &result.termCandidates);
        collectArrays(first->GetCellData(), ArrayInfo::Association::Cell,
                      &arrays, &result.termCandidates);
        commonFields = actualArrayNames(first);
        if (commonFields.isEmpty()) {
            result.errorMessage = "The first EnSight Gold frame does not contain numeric point or cell arrays.";
            return result;
        }
    }

    if (validateEveryFrame) {
        for (int index = 1; index < frames.size(); ++index) {
            if (canceled()) {
                return result;
            }
            const EnSightGoldFrameSource& frameSource = frames.at(index);
            QString frameError;
            vtkSmartPointer<vtkDataSet> frame = readEnSightGoldVolumeFromCaseAtTime(cases.at(frameSource.caseIndex).path,
                                                                                    frameSource.timeValue,
                                                                                    &frameError);
            if (frame == nullptr) {
                result.errorMessage = frameError;
                return result;
            }
            if (coordinateCheck && coordinateSignature(frame) != firstCoordinateSignature &&
                result.warningMessage.isEmpty()) {
                result.warningMessage = QString("EnSight Gold frame coordinates differ from the first frame at time %1 in %2.\nFirst bounds: %3\nCurrent bounds: %4\nImport will continue; verify the combined dataset before analysis.")
                                            .arg(formatDouble(frameSource.timeValue),
                                                 cases.at(frameSource.caseIndex).path,
                                                 firstBoundsText,
                                                 boundsSummary(frame));
            }
            if (commonFieldsOnly) {
                commonFields = intersectFieldNames(commonFields, actualArrayNames(frame));
                if (commonFields.isEmpty()) {
                    result.errorMessage = QStringLiteral("No numeric EnSight Gold variables are present in every time step.");
                    return result;
                }
            }
            if (progress) {
                progress(QStringLiteral("Validating EnSight Gold case"), index + 1, frames.size());
            }
            if (canceled()) {
                return result;
            }
        }
    } else if (progress) {
        progress(QStringLiteral("Validating EnSight Gold case"), 2, 2);
        if (canceled()) {
            return result;
        }
    }

    if (commonFieldsOnly) {
        arrays = filterArrayInfos(arrays, commonFields);
    }
    result.termCandidates = candidateListForArrays(arrays);
    inferTerms(arrays, dimension, &result.entry);
    result.entry.termCandidates = result.termCandidates;

    if (result.termCandidates.isEmpty()) {
        result.errorMessage = commonFieldsOnly
            ? QStringLiteral("No numeric EnSight Gold variables are present in every time step.")
            : QStringLiteral("The first EnSight Gold frame does not contain numeric point or cell arrays.");
        return result;
    }

    result.entry.checkStatus = QStringLiteral("passed");
    result.passed = true;
    return result;
}

bool FlowDataInspector::convertOpenFoamCaseToPvd(const QString& casePath,
                                                 const QString& caseType,
                                                 const QString& destinationFolder,
                                                 const QString& pvdBaseName,
                                                 bool startEndFieldsOnly,
                                                 const QStringList& allowedFields,
                                                 QString* errorMessage,
                                                 const ProgressCallback& progress,
                                                 const CancellationCallback& shouldCancel) {
    auto fail = [errorMessage](const QString& message) {
        if (errorMessage != nullptr) {
            *errorMessage = message;
        }
        return false;
    };
    auto canceled = [&]() {
        return shouldCancel && shouldCancel();
    };

    if (destinationFolder.trimmed().isEmpty()) {
        return fail(QStringLiteral("OpenFOAM conversion destination is empty."));
    }
    const QStringList allowedActualFields = actualFieldsFromCandidates(allowedFields);
    if (startEndFieldsOnly && allowedActualFields.isEmpty()) {
        return fail(QStringLiteral("OpenFOAM conversion requires validated common fields."));
    }
    if (canceled()) {
        return fail(QStringLiteral("OpenFOAM conversion interrupted."));
    }

    vtkSmartPointer<vtkOpenFOAMReader> reader = createOpenFoamReader(caseType);
    QString readerError;
    if (!prepareOpenFoamReader(reader, casePath, caseType, &readerError)) {
        return fail(readerError);
    }

    const QVector<double> timeValues = openFoamTimeValues(reader);
    if (timeValues.isEmpty()) {
        return fail(QStringLiteral("The OpenFOAM case does not contain any time steps."));
    }
    if (!QDir().mkpath(destinationFolder)) {
        return fail(QString("Cannot create OpenFOAM conversion destination: %1").arg(destinationFolder));
    }

    const QString baseName = pvdBaseName.trimmed().isEmpty() ? QStringLiteral("openfoam") : pvdBaseName.trimmed();
    QStringList relativeFramePaths;
    relativeFramePaths.reserve(timeValues.size());
    for (int index = 0; index < timeValues.size(); ++index) {
        if (canceled()) {
            return fail(QStringLiteral("OpenFOAM conversion interrupted."));
        }
        QString frameError;
        vtkSmartPointer<vtkDataSet> frame = readOpenFoamVolumeAtTime(reader, timeValues.at(index), &frameError);
        if (frame == nullptr) {
            return fail(frameError);
        }
        if (startEndFieldsOnly) {
            filterDataSetArrays(frame, allowedActualFields);
        }

        const QString frameName = QStringLiteral("%1_%2.vtu")
                                      .arg(baseName)
                                      .arg(index, 6, 10, QLatin1Char('0'));
        const QString framePath = normalizedPath(QDir(destinationFolder).absoluteFilePath(frameName));
        if (!writeUnstructuredFrame(frame, framePath, &frameError)) {
            return fail(frameError);
        }
        relativeFramePaths.append(frameName);
        if (progress) {
            progress(QStringLiteral("Converting OpenFOAM case"), index + 1, timeValues.size() + 1);
        }
    }

    const QString pvdPath = normalizedPath(QDir(destinationFolder).absoluteFilePath(baseName + QStringLiteral(".pvd")));
    QString pvdError;
    if (!writePvdFile(pvdPath, timeValues, relativeFramePaths, &pvdError)) {
        return fail(pvdError);
    }
    if (progress) {
        progress(QStringLiteral("Converting OpenFOAM case"), timeValues.size() + 1, timeValues.size() + 1);
    }
    return true;
}

bool FlowDataInspector::convertEnSightGoldCaseToPvd(const QString& casePath,
                                                    const QString& destinationFolder,
                                                    const QString& pvdBaseName,
                                                    bool commonFieldsOnly,
                                                    const QStringList& allowedFields,
                                                    QString* errorMessage,
                                                    const ProgressCallback& progress,
                                                    const CancellationCallback& shouldCancel) {
    return convertEnSightGoldCasesToPvd({casePath},
                                        destinationFolder,
                                        pvdBaseName,
                                        commonFieldsOnly,
                                        allowedFields,
                                        errorMessage,
                                        progress,
                                        shouldCancel);
}

bool FlowDataInspector::convertEnSightGoldCasesToPvd(const QStringList& casePaths,
                                                     const QString& destinationFolder,
                                                     const QString& pvdBaseName,
                                                     bool commonFieldsOnly,
                                                     const QStringList& allowedFields,
                                                     QString* errorMessage,
                                                     const ProgressCallback& progress,
                                                     const CancellationCallback& shouldCancel) {
    auto fail = [errorMessage](const QString& message) {
        if (errorMessage != nullptr) {
            *errorMessage = message;
        }
        return false;
    };
    auto canceled = [&]() {
        return shouldCancel && shouldCancel();
    };

    if (destinationFolder.trimmed().isEmpty()) {
        return fail(QStringLiteral("EnSight Gold conversion destination is empty."));
    }
    const QStringList allowedActualFields = actualFieldsFromCandidates(allowedFields);
    if (commonFieldsOnly && allowedActualFields.isEmpty()) {
        return fail(QStringLiteral("EnSight Gold conversion requires validated common variables."));
    }
    if (canceled()) {
        return fail(QStringLiteral("EnSight Gold conversion interrupted."));
    }

    QVector<EnSightGoldCaseSource> cases;
    QVector<EnSightGoldFrameSource> frames;
    QString readerError;
    if (!loadEnSightGoldCaseSources(casePaths, &cases, &frames, &readerError)) {
        return fail(readerError);
    }
    if (frames.isEmpty()) {
        return fail(QStringLiteral("The EnSight Gold case does not contain any readable frames."));
    }
    QString geometryCompatibilityError;
    if (!checkEnSightGoldSourcesForVtk(cases, &geometryCompatibilityError)) {
        return fail(geometryCompatibilityError);
    }
    if (!QDir().mkpath(destinationFolder)) {
        return fail(QString("Cannot create EnSight Gold conversion destination: %1").arg(destinationFolder));
    }

    const QString baseName = pvdBaseName.trimmed().isEmpty() ? QStringLiteral("ensight_gold") : pvdBaseName.trimmed();
    QVector<double> timeValues;
    timeValues.reserve(frames.size());
    QStringList relativeFramePaths;
    relativeFramePaths.reserve(frames.size());
    for (int index = 0; index < frames.size(); ++index) {
        if (canceled()) {
            return fail(QStringLiteral("EnSight Gold conversion interrupted."));
        }
        const EnSightGoldFrameSource& frameSource = frames.at(index);
        QString frameError;
        vtkSmartPointer<vtkDataSet> frame;
        try {
            frame = readEnSightGoldVolumeFromCaseAtTime(cases.at(frameSource.caseIndex).path,
                                                        frameSource.timeValue,
                                                        &frameError);
        } catch (const std::exception& ex) {
            return fail(QString("VTK failed while reading EnSight Gold frame %1 at time %2 from %3: %4")
                            .arg(index + 1)
                            .arg(formatDouble(frameSource.timeValue),
                                 cases.at(frameSource.caseIndex).path,
                                 QString::fromLocal8Bit(ex.what())));
        } catch (...) {
            return fail(QString("VTK failed with an unknown exception while reading EnSight Gold frame %1 at time %2 from %3.")
                            .arg(index + 1)
                            .arg(formatDouble(frameSource.timeValue),
                                 cases.at(frameSource.caseIndex).path));
        }
        if (frame == nullptr) {
            return fail(frameError);
        }
        if (commonFieldsOnly) {
            filterDataSetArrays(frame, allowedActualFields);
        }

        const QString frameName = QStringLiteral("%1_%2.vtu")
                                      .arg(baseName)
                                      .arg(index, 6, 10, QLatin1Char('0'));
        const QString framePath = normalizedPath(QDir(destinationFolder).absoluteFilePath(frameName));
        try {
            if (!writeUnstructuredFrame(frame, framePath, &frameError)) {
                return fail(frameError);
            }
        } catch (const std::exception& ex) {
            return fail(QString("VTK failed while writing converted EnSight Gold frame %1 to %2: %3")
                            .arg(index + 1)
                            .arg(framePath,
                                 QString::fromLocal8Bit(ex.what())));
        } catch (...) {
            return fail(QString("VTK failed with an unknown exception while writing converted EnSight Gold frame %1 to %2.")
                            .arg(index + 1)
                            .arg(framePath));
        }
        timeValues.append(frameSource.timeValue);
        relativeFramePaths.append(frameName);
        if (progress) {
            progress(QStringLiteral("Converting EnSight Gold case"), index + 1, frames.size() + 1);
        }
    }

    const QString pvdPath = normalizedPath(QDir(destinationFolder).absoluteFilePath(baseName + QStringLiteral(".pvd")));
    QString pvdError;
    if (!writePvdFile(pvdPath, timeValues, relativeFramePaths, &pvdError)) {
        return fail(pvdError);
    }
    if (progress) {
        progress(QStringLiteral("Converting EnSight Gold case"), frames.size() + 1, frames.size() + 1);
    }
    return true;
}
