#include "FlowDataInspector.h"
#include "streamcenterplus_resample/Core.h"

#include <streamcenterplus/VtkHdfIO.h>

#include <vtkAbstractArray.h>
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
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
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QRegularExpression>
#include <QSaveFile>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

namespace {

QString normalizedPath(const QString& path) {
    return QDir::fromNativeSeparators(path);
}

std::filesystem::path fsPathFromQString(const QString& path) {
    return std::filesystem::path(QDir::toNativeSeparators(path).toStdWString());
}

QString qStringFromFsPath(const std::filesystem::path& path) {
    return normalizedPath(QString::fromStdWString(path.wstring()));
}

QString realVtkPath(const QString& path) {
    return qStringFromFsPath(streamcenterplus::vtkhdf::RealPath(fsPathFromQString(path)));
}

QString vtkExtensionSuffix(const QString& path) {
    std::string ext = streamcenterplus::vtkhdf::Extension(fsPathFromQString(path));
    if (!ext.empty() && ext.front() == '.') {
        ext.erase(ext.begin());
    }
    return QString::fromStdString(ext);
}

QString formatDouble(double value) {
    return QString::number(value, 'g', 15);
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
        || (includeNek5000 &&
            (ext == QStringLiteral("nek5000") || ext == QStringLiteral("nek")));
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

QStringList collectVtkFiles(const QString& folder, bool includeNek5000) {
    QDir dir(folder);
    const QFileInfoList entries = dir.entryInfoList(QDir::Files | QDir::NoDotAndDotDot | QDir::Readable,
                                                    QDir::Name | QDir::IgnoreCase);
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
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1e-12;
                 }),
                 values.end());
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
    if (!caseInfo.exists() || !caseInfo.isFile()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("EnSight Gold case file does not exist: %1").arg(casePath);
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
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1e-12;
                 }),
                 values.end());
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
        xml.writeAttribute(QStringLiteral("timestep"), formatDouble(timeValues.at(index)));
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

bool validateNek5000Descriptor(const QString& path, QString* errorMessage) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot open Nek5000 descriptor: %1").arg(file.errorString());
        }
        return false;
    }

    bool hasFileTemplate = false;
    bool hasFirstTimeStep = false;
    bool hasNumberOfTimeSteps = false;
    const QStringList lines = QString::fromUtf8(file.readAll()).split(QLatin1Char('\n'));
    for (QString line : lines) {
        line = line.trimmed();
        if (line.isEmpty() || line.startsWith(QLatin1Char('#'))) {
            continue;
        }

        int tokenEnd = 0;
        while (tokenEnd < line.size() && !line.at(tokenEnd).isSpace()) {
            ++tokenEnd;
        }
        const QString token = line.left(tokenEnd).toLower();
        if (token == QStringLiteral("nek5000")
            || token == QStringLiteral("endian:")
            || token == QStringLiteral("version:")) {
            continue;
        }
        if (token == QStringLiteral("filetemplate:")) {
            hasFileTemplate = true;
            continue;
        }
        if (token == QStringLiteral("firsttimestep:")) {
            hasFirstTimeStep = true;
            continue;
        }
        if (token == QStringLiteral("numtimesteps:")) {
            hasNumberOfTimeSteps = true;
            continue;
        }

        if (errorMessage != nullptr) {
            *errorMessage = QString("Unsupported Nek5000 descriptor tag \"%1\" in %2.").arg(token, path);
        }
        return false;
    }

    if (!hasFileTemplate || !hasFirstTimeStep || !hasNumberOfTimeSteps) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Nek5000 descriptor must define filetemplate:, firsttimestep:, and numtimesteps:: %1").arg(path);
        }
        return false;
    }
    return true;
}

QVector<double> nek5000TimeValues(const QString& filePath, QString* errorMessage) {
    if (!validateNek5000Descriptor(filePath, errorMessage)) {
        return {};
    }

    QVector<double> values;
    try {
        const std::vector<streamcenterplus::resample::FrameInfo> frames =
            streamcenterplus::resample::EnumerateDataSetFrames(fsPathFromQString(filePath));
        values.reserve(static_cast<int>(std::min<std::size_t>(
            frames.size(), static_cast<std::size_t>(std::numeric_limits<int>::max()))));
        for (const streamcenterplus::resample::FrameInfo& frame : frames) {
            values.push_back(frame.time);
        }
    } catch (const std::exception& ex) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot enumerate Nek5000 time steps from %1:\n%2")
                                .arg(filePath, QString::fromLocal8Bit(ex.what()));
        }
        return {};
    }

    std::sort(values.begin(), values.end());
    QVector<double> uniqueValues;
    uniqueValues.reserve(values.size());
    for (const double value : std::as_const(values)) {
        if (uniqueValues.isEmpty() || std::abs(uniqueValues.back() - value) > 1.0e-10) {
            uniqueValues.push_back(value);
        }
    }
    return uniqueValues;
}

vtkSmartPointer<vtkDataSet> readAnyVtk(const QString& filePath, QString* errorMessage) {
    const QString ext = vtkExtensionSuffix(filePath);
    vtkDataSet* dataSet = nullptr;

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

    const QString realPath = realVtkPath(filePath);
    const QByteArray realPathBytes = QFile::encodeName(realPath);

    if (ext == QStringLiteral("nek5000") || ext == QStringLiteral("nek")) {
        if (!validateNek5000Descriptor(filePath, errorMessage)) {
            return nullptr;
        }
        try {
            const std::vector<streamcenterplus::resample::FrameInfo> frames =
                streamcenterplus::resample::EnumerateDataSetFrames(fsPathFromQString(filePath));
            if (!frames.empty()) {
                return streamcenterplus::resample::ReadDataSet(frames.front());
            }
        } catch (const std::exception& ex) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("Cannot read Nek5000 dataset %1:\n%2")
                                    .arg(filePath, QString::fromLocal8Bit(ex.what()));
            }
            return nullptr;
        }
    } else if (ext == QStringLiteral("vtk")) {
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
    if (values.size() <= 2) {
        return true;
    }
    const double spacing = values[1] - values[0];
    const double tolerance = std::max(1.0, std::abs(spacing)) * 1e-7;
    for (int index = 2; index < values.size(); ++index) {
        if (std::abs((values[index] - values[index - 1]) - spacing) > tolerance) {
            return false;
        }
    }
    return true;
}

QVector<double> uniqueAxisValues(vtkDataSet* dataSet, int axis) {
    QVector<double> values;
    values.reserve(static_cast<int>(std::min<vtkIdType>(dataSet->GetNumberOfPoints(), std::numeric_limits<int>::max())));
    double point[3] = {};
    for (vtkIdType id = 0; id < dataSet->GetNumberOfPoints(); ++id) {
        dataSet->GetPoint(id, point);
        values.append(point[axis]);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1e-10;
                 }),
                 values.end());
    return values;
}

bool isUniformStructured(vtkDataSet* dataSet, int dims[3]) {
    if (vtkImageData::SafeDownCast(dataSet) != nullptr) {
        return dimensionsForDataSet(dataSet, dims);
    }
    if (!dimensionsForDataSet(dataSet, dims)) {
        return false;
    }
    return axisIsRegular(uniqueAxisValues(dataSet, 0))
        && axisIsRegular(uniqueAxisValues(dataSet, 1))
        && axisIsRegular(uniqueAxisValues(dataSet, 2));
}

struct ArrayInfo {
    QString name;
    int components = 0;
    vtkIdType tuples = 0;
};

void collectArrays(vtkFieldData* data, QList<ArrayInfo>* arrays, QStringList* candidates) {
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

struct VtkDataSummary {
    bool hasBounds = false;
    double bounds[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int dims[3] = {0, 0, 0};
    bool hasDims = false;
    bool uniformStructured = false;
    QList<ArrayInfo> arrays;
    QStringList termCandidates;
};

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
    if (dataSet == nullptr || summary == nullptr
        || dataSet->GetNumberOfPoints() <= 0
        || dataSet->GetNumberOfCells() <= 0) {
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
    summary->uniformStructured = summary->uniformStructured || isUniformStructured(dataSet, dims);

    if (summary->arrays.isEmpty()) {
        collectArrays(dataSet->GetPointData(), &summary->arrays, &summary->termCandidates);
        collectArrays(dataSet->GetCellData(), &summary->arrays, &summary->termCandidates);
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

bool boundsEquivalent(const VtkDataSummary& lhs, const VtkDataSummary& rhs) {
    if (!lhs.hasBounds || !rhs.hasBounds) {
        return false;
    }
    for (int index = 0; index < 6; ++index) {
        const double scale = std::max({1.0, std::abs(lhs.bounds[index]), std::abs(rhs.bounds[index])});
        if (std::abs(lhs.bounds[index] - rhs.bounds[index]) > scale * 1e-8) {
            return false;
        }
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
    return true;
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

bool boundsEquivalent(vtkDataSet* lhs, vtkDataSet* rhs) {
    double left[6] = {};
    double right[6] = {};
    lhs->GetBounds(left);
    rhs->GetBounds(right);
    for (int index = 0; index < 6; ++index) {
        const double scale = std::max({1.0, std::abs(left[index]), std::abs(right[index])});
        if (std::abs(left[index] - right[index]) > scale * 1e-8) {
            return false;
        }
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
    return true;
}

bool boundsRangeEquivalent(const double* left, vtkDataSet* rhs) {
    double right[6] = {};
    rhs->GetBounds(right);
    for (int index = 0; index < 6; ++index) {
        const double scale = std::max({1.0, std::abs(left[index]), std::abs(right[index])});
        if (std::abs(left[index] - right[index]) > scale * 1e-8) {
            return false;
        }
    }
    return true;
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

QStringList splitEnSightGoldSources(const QStringList& sources) {
    QStringList expanded;
    for (const QString& source : sources) {
        const QString trimmed = source.trimmed();
        if (trimmed.isEmpty()) {
            continue;
        }
        if (QFileInfo(trimmed).exists()) {
            expanded.append(trimmed);
            continue;
        }
        const QStringList pieces = trimmed.split(QLatin1Char(';'), Qt::SkipEmptyParts);
        for (const QString& piece : pieces) {
            const QString path = piece.trimmed();
            if (!path.isEmpty()) {
                expanded.append(path);
            }
        }
    }
    return expanded;
}

QStringList collectEnSightGoldCaseFiles(const QStringList& sources, QString* errorMessage) {
    QStringList caseFiles;
    const QStringList expandedSources = splitEnSightGoldSources(sources);
    if (expandedSources.isEmpty()) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Select at least one EnSight Gold *.case file or a folder containing *.case files.");
        }
        return caseFiles;
    }

    for (const QString& source : expandedSources) {
        const QFileInfo info(source);
        if (!info.exists()) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold source does not exist: %1").arg(source);
            }
            caseFiles.clear();
            return caseFiles;
        }

        if (info.isDir()) {
            const QFileInfoList entries = QDir(info.absoluteFilePath()).entryInfoList({QStringLiteral("*.case")},
                                                                                     QDir::Files | QDir::NoDotAndDotDot | QDir::Readable,
                                                                                     QDir::Name | QDir::IgnoreCase);
            if (entries.isEmpty()) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("No EnSight Gold case files (*.case) were found in folder: %1")
                                        .arg(info.absoluteFilePath());
                }
                caseFiles.clear();
                return caseFiles;
            }
            for (const QFileInfo& entry : entries) {
                const QString path = normalizedPath(entry.absoluteFilePath());
                if (!caseFiles.contains(path, Qt::CaseInsensitive)) {
                    caseFiles.append(path);
                }
            }
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
        if (!caseFiles.contains(path, Qt::CaseInsensitive)) {
            caseFiles.append(path);
        }
    }

    return caseFiles;
}

struct EnSightGoldCaseSource {
    QString path;
    QVector<double> timeValues;
    QList<ArrayInfo> arrays;
    QStringList geometryPatterns;
    QStringList variablePatterns;
    int filenameStartNumber = 0;
    int filenameIncrement = 1;
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

bool tokenIsInteger(const QString& text) {
    int ignored = 0;
    return parseIntToken(text, &ignored);
}

QVector<double> parseDoubleTokens(QString text) {
    QVector<double> values;
    text.replace(QLatin1Char(','), QLatin1Char(' '));
    const QStringList tokens = text.split(QRegularExpression(QStringLiteral(R"(\s+)")), Qt::SkipEmptyParts);
    values.reserve(tokens.size());
    for (const QString& token : tokens) {
        bool ok = false;
        const double value = token.toDouble(&ok);
        if (ok && std::isfinite(value)) {
            values.push_back(value);
        }
    }
    return values;
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

QStringList enSightGoldPayloadTokens(const QString& payload) {
    return payload.trimmed().split(QRegularExpression(QStringLiteral(R"(\s+)")), Qt::SkipEmptyParts);
}

bool appendEnSightGoldGeometryPattern(const QString& payload, EnSightGoldCaseSource* source) {
    QStringList tokens = enSightGoldPayloadTokens(payload);
    if (tokens.isEmpty() || source == nullptr) {
        return false;
    }
    if (tokens.size() >= 2 && tokenIsInteger(tokens.first())) {
        tokens.removeFirst();
    }
    const QString pattern = cleanEnSightGoldName(tokens.join(QLatin1Char(' ')));
    if (!pattern.isEmpty() && !source->geometryPatterns.contains(pattern, Qt::CaseInsensitive)) {
        source->geometryPatterns.append(pattern);
    }
    return !pattern.isEmpty();
}

bool appendEnSightGoldVariable(const QString& declaration,
                               const QString& payload,
                               EnSightGoldCaseSource* source) {
    QStringList tokens = enSightGoldPayloadTokens(payload);
    if (tokens.size() < 2 || source == nullptr) {
        return false;
    }
    if (tokens.size() >= 3 && tokenIsInteger(tokens.first())) {
        tokens.removeFirst();
    }

    const QString pattern = cleanEnSightGoldName(tokens.takeLast());
    const QString name = cleanEnSightGoldName(tokens.join(QLatin1Char(' ')));
    if (name.isEmpty() || pattern.isEmpty()) {
        return false;
    }

    bool alreadyPresent = false;
    for (const ArrayInfo& array : std::as_const(source->arrays)) {
        if (array.name.compare(name, Qt::CaseInsensitive) == 0) {
            alreadyPresent = true;
            break;
        }
    }
    if (!alreadyPresent) {
        ArrayInfo info;
        info.name = name;
        info.components = enSightGoldComponentCount(declaration);
        source->arrays.append(info);
    }
    if (!source->variablePatterns.contains(pattern, Qt::CaseInsensitive)) {
        source->variablePatterns.append(pattern);
    }
    return true;
}

QString enSightGoldSectionName(const QString& line) {
    const QString upper = line.trimmed().toUpper();
    if (upper == QStringLiteral("FORMAT")
        || upper == QStringLiteral("GEOMETRY")
        || upper == QStringLiteral("VARIABLE")
        || upper == QStringLiteral("TIME")) {
        return upper;
    }
    return {};
}

bool parseEnSightGoldCaseSource(const QString& casePath,
                                EnSightGoldCaseSource* source,
                                QString* errorMessage) {
    if (source == nullptr) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("Internal error: EnSight Gold metadata storage was not provided.");
        }
        return false;
    }

    const QFileInfo caseInfo(casePath);
    if (!caseInfo.exists() || !caseInfo.isFile()) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("EnSight Gold case file does not exist: %1").arg(casePath);
        }
        return false;
    }
    if (caseInfo.suffix().compare(QStringLiteral("case"), Qt::CaseInsensitive) != 0) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("EnSight Gold import expects a *.case file.");
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
    source->path = normalizedPath(caseInfo.absoluteFilePath());

    const QString text = QString::fromUtf8(file.readAll());
    QString section;
    bool foundGoldFormat = false;
    bool readingTimeValues = false;
    int declaredSteps = 0;
    const QStringList lines = text.split(QRegularExpression(QStringLiteral(R"(\r\n|\n|\r)")));
    for (const QString& rawLine : lines) {
        QString line = rawLine.trimmed();
        if (line.isEmpty() || line.startsWith(QLatin1Char('#'))) {
            continue;
        }

        const QString nextSection = enSightGoldSectionName(line);
        if (!nextSection.isEmpty()) {
            section = nextSection;
            readingTimeValues = false;
            continue;
        }

        const int colon = line.indexOf(QLatin1Char(':'));
        if (section == QStringLiteral("TIME") && readingTimeValues && colon < 0) {
            const QVector<double> values = parseDoubleTokens(line);
            source->timeValues += values;
            continue;
        }
        if (colon < 0) {
            continue;
        }

        const QString key = line.left(colon).trimmed().toLower();
        const QString value = line.mid(colon + 1).trimmed();
        readingTimeValues = false;

        if (section == QStringLiteral("FORMAT")) {
            if (key == QStringLiteral("type") && value.toLower().contains(QStringLiteral("ensight gold"))) {
                foundGoldFormat = true;
            }
            continue;
        }
        if (section == QStringLiteral("GEOMETRY")) {
            appendEnSightGoldGeometryPattern(value, source);
            continue;
        }
        if (section == QStringLiteral("VARIABLE")) {
            appendEnSightGoldVariable(key, value, source);
            continue;
        }
        if (section == QStringLiteral("TIME")) {
            if (key == QStringLiteral("number of steps")) {
                parseIntToken(value, &declaredSteps);
            } else if (key == QStringLiteral("filename start number")) {
                parseIntToken(value, &source->filenameStartNumber);
            } else if (key == QStringLiteral("filename increment")) {
                int increment = 0;
                if (parseIntToken(value, &increment) && increment != 0) {
                    source->filenameIncrement = increment;
                }
            } else if (key == QStringLiteral("time values")) {
                source->timeValues += parseDoubleTokens(value);
                readingTimeValues = true;
            }
        }
    }

    if (!foundGoldFormat) {
        if (errorMessage != nullptr) {
            *errorMessage = QStringLiteral("The selected *.case file is not an EnSight Gold case.");
        }
        return false;
    }
    if (source->geometryPatterns.isEmpty()) {
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
    if (source->timeValues.isEmpty() && declaredSteps > 0) {
        source->timeValues.reserve(declaredSteps);
        for (int index = 0; index < declaredSteps; ++index) {
            source->timeValues.push_back(static_cast<double>(index));
        }
    }
    if (source->timeValues.isEmpty()) {
        source->timeValues.push_back(0.0);
    }
    return true;
}

QString expandEnSightGoldPattern(QString pattern, int frameIndex, const EnSightGoldCaseSource& source) {
    const QRegularExpression stars(QStringLiteral(R"(\*+)"));
    const QRegularExpressionMatch match = stars.match(pattern);
    if (!match.hasMatch()) {
        return pattern;
    }

    const int width = match.capturedLength();
    const int frameNumber = source.filenameStartNumber + frameIndex * source.filenameIncrement;
    const QString number = QString::number(frameNumber).rightJustified(width, QLatin1Char('0'));
    pattern.replace(match.capturedStart(), width, number);
    return pattern;
}

QString resolveEnSightGoldReferencedPath(const EnSightGoldCaseSource& source,
                                         const QString& pattern,
                                         int frameIndex) {
    const QString expanded = expandEnSightGoldPattern(pattern, frameIndex, source);
    if (QDir::isAbsolutePath(expanded)) {
        return normalizedPath(QFileInfo(expanded).absoluteFilePath());
    }
    return normalizedPath(QFileInfo(QDir(QFileInfo(source.path).absolutePath()).absoluteFilePath(expanded)).absoluteFilePath());
}

bool checkEnSightGoldReference(const EnSightGoldCaseSource& source,
                               const QString& pattern,
                               int frameIndex,
                               QString* errorMessage) {
    const QString path = resolveEnSightGoldReferencedPath(source, pattern, frameIndex);
    if (QFileInfo::exists(path)) {
        return true;
    }
    if (errorMessage != nullptr) {
        *errorMessage = QString("EnSight Gold case references a missing file: %1").arg(path);
    }
    return false;
}

bool verifyEnSightGoldReferences(const EnSightGoldCaseSource& source, QString* errorMessage) {
    const int lastFrameIndex = source.timeValues.isEmpty() ? 0 : static_cast<int>(source.timeValues.size() - 1);
    auto checkFirstAndLast = [&](const QString& pattern) {
        if (!checkEnSightGoldReference(source, pattern, 0, errorMessage)) {
            return false;
        }
        if (lastFrameIndex > 0 && pattern.contains(QLatin1Char('*'))) {
            return checkEnSightGoldReference(source, pattern, lastFrameIndex, errorMessage);
        }
        return true;
    };

    for (const QString& pattern : source.geometryPatterns) {
        if (!checkFirstAndLast(pattern)) {
            return false;
        }
    }
    for (const QString& pattern : source.variablePatterns) {
        if (!checkFirstAndLast(pattern)) {
            return false;
        }
    }
    return true;
}

bool enSightGoldSourceLooksLarge(const EnSightGoldCaseSource& source) {
    QStringList firstFramePaths;
    firstFramePaths.reserve(source.geometryPatterns.size() + source.variablePatterns.size());
    for (const QString& pattern : source.geometryPatterns) {
        const QString path = resolveEnSightGoldReferencedPath(source, pattern, 0);
        if (!firstFramePaths.contains(path, Qt::CaseInsensitive)) {
            firstFramePaths.append(path);
        }
    }
    for (const QString& pattern : source.variablePatterns) {
        const QString path = resolveEnSightGoldReferencedPath(source, pattern, 0);
        if (!firstFramePaths.contains(path, Qt::CaseInsensitive)) {
            firstFramePaths.append(path);
        }
    }

    qint64 firstFrameBytes = 0;
    for (const QString& path : std::as_const(firstFramePaths)) {
        const qint64 bytes = QFileInfo(path).size();
        if (bytes >= kLargeEnSightGoldSingleFileBytes) {
            return true;
        }
        if (bytes > 0) {
            firstFrameBytes += bytes;
        }
        if (firstFrameBytes >= kLargeEnSightGoldFrameBytes) {
            return true;
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

bool readEnSightGoldTextRecord(QFile* file, QString* text) {
    if (file == nullptr || text == nullptr) {
        return false;
    }
    const QByteArray bytes = file->read(80);
    if (bytes.size() != 80) {
        return false;
    }
    *text = QString::fromLatin1(bytes).trimmed();
    return true;
}

bool readEnSightGoldInt32(QFile* file, qint32* value) {
    if (file == nullptr || value == nullptr) {
        return false;
    }
    const QByteArray bytes = file->read(4);
    if (bytes.size() != 4) {
        return false;
    }
    const uchar* raw = reinterpret_cast<const uchar*>(bytes.constData());
    *value = static_cast<qint32>(raw[0])
        | (static_cast<qint32>(raw[1]) << 8)
        | (static_cast<qint32>(raw[2]) << 16)
        | (static_cast<qint32>(raw[3]) << 24);
    return true;
}

bool skipEnSightGoldBytes(QFile* file, qint64 byteCount) {
    if (file == nullptr || byteCount < 0) {
        return false;
    }
    return file->seek(file->pos() + byteCount);
}

int enSightGoldFixedNodesPerElement(const QString& elementType) {
    const QString type = elementType.trimmed().toLower();
    if (type == QStringLiteral("point")) {
        return 1;
    }
    if (type == QStringLiteral("bar2")) {
        return 2;
    }
    if (type == QStringLiteral("bar3") || type == QStringLiteral("tria3")) {
        return 3;
    }
    if (type == QStringLiteral("quad4") || type == QStringLiteral("tetra4")) {
        return 4;
    }
    if (type == QStringLiteral("pyramid5")) {
        return 5;
    }
    if (type == QStringLiteral("tria6") || type == QStringLiteral("penta6")) {
        return 6;
    }
    if (type == QStringLiteral("quad8") || type == QStringLiteral("hexa8")) {
        return 8;
    }
    if (type == QStringLiteral("tetra10")) {
        return 10;
    }
    if (type == QStringLiteral("pyramid13")) {
        return 13;
    }
    if (type == QStringLiteral("penta15")) {
        return 15;
    }
    if (type == QStringLiteral("hexa20")) {
        return 20;
    }
    return 0;
}

bool sumEnSightGoldInt32Block(QFile* file,
                              qint64 count,
                              qint64* sum,
                              qint32* minValue,
                              qint32* maxValue) {
    if (file == nullptr || count < 0 || sum == nullptr) {
        return false;
    }
    constexpr qint64 kChunkInts = 1024 * 1024;
    *sum = 0;
    if (minValue != nullptr) {
        *minValue = std::numeric_limits<qint32>::max();
    }
    if (maxValue != nullptr) {
        *maxValue = std::numeric_limits<qint32>::min();
    }

    QByteArray buffer;
    buffer.resize(static_cast<int>(kChunkInts * 4));
    qint64 consumed = 0;
    while (consumed < count) {
        const qint64 intsToRead = std::min(kChunkInts, count - consumed);
        const qint64 bytesToRead = intsToRead * 4;
        const qint64 bytesRead = file->read(buffer.data(), bytesToRead);
        if (bytesRead != bytesToRead) {
            return false;
        }
        const uchar* raw = reinterpret_cast<const uchar*>(buffer.constData());
        for (qint64 index = 0; index < intsToRead; ++index) {
            const qint64 offset = index * 4;
            const qint32 value = static_cast<qint32>(raw[offset])
                | (static_cast<qint32>(raw[offset + 1]) << 8)
                | (static_cast<qint32>(raw[offset + 2]) << 16)
                | (static_cast<qint32>(raw[offset + 3]) << 24);
            if (value < 0) {
                return false;
            }
            *sum += value;
            if (minValue != nullptr) {
                *minValue = std::min(*minValue, value);
            }
            if (maxValue != nullptr) {
                *maxValue = std::max(*maxValue, value);
            }
        }
        consumed += intsToRead;
    }
    return true;
}

bool checkEnSightGoldGeometryForVtk(const QString& geometryPath, QString* errorMessage) {
    QFile file(geometryPath);
    if (!file.open(QIODevice::ReadOnly)) {
        if (errorMessage != nullptr) {
            *errorMessage = QString("Cannot inspect EnSight Gold geometry file: %1").arg(geometryPath);
        }
        return false;
    }

    QString text;
    for (int header = 0; header < 5; ++header) {
        if (!readEnSightGoldTextRecord(&file, &text)) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold geometry header is incomplete: %1").arg(geometryPath);
            }
            return false;
        }
    }

    while (file.pos() < file.size()) {
        if (!readEnSightGoldTextRecord(&file, &text)) {
            break;
        }
        if (text != QStringLiteral("part")) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("Unexpected EnSight Gold geometry record \"%1\" in %2 at byte %3.")
                                    .arg(text, geometryPath)
                                    .arg(file.pos() - 80);
            }
            return false;
        }

        qint32 partId = 0;
        qint32 nodeCount = 0;
        if (!readEnSightGoldInt32(&file, &partId)
            || !readEnSightGoldTextRecord(&file, &text)
            || !readEnSightGoldTextRecord(&file, &text)
            || text != QStringLiteral("coordinates")
            || !readEnSightGoldInt32(&file, &nodeCount)
            || nodeCount < 0) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold geometry part metadata is invalid in %1.").arg(geometryPath);
            }
            return false;
        }
        if (!skipEnSightGoldBytes(&file, static_cast<qint64>(nodeCount) * 3 * 4)) {
            if (errorMessage != nullptr) {
                *errorMessage = QString("EnSight Gold geometry coordinates are truncated in %1, part %2.")
                                    .arg(geometryPath)
                                    .arg(partId);
            }
            return false;
        }

        while (file.pos() < file.size()) {
            const qint64 recordOffset = file.pos();
            if (!readEnSightGoldTextRecord(&file, &text)) {
                break;
            }
            if (text == QStringLiteral("part")) {
                file.seek(recordOffset);
                break;
            }

            qint32 elementCount = 0;
            if (!readEnSightGoldInt32(&file, &elementCount) || elementCount < 0) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("EnSight Gold geometry element count is invalid in %1, part %2, element type %3.")
                                        .arg(geometryPath)
                                        .arg(partId)
                                        .arg(text);
                }
                return false;
            }

            if (text == QStringLiteral("nfaced")) {
                qint64 faceCount = 0;
                qint32 minFaces = 0;
                qint32 maxFaces = 0;
                if (!sumEnSightGoldInt32Block(&file, elementCount, &faceCount, &minFaces, &maxFaces)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QString("EnSight Gold nfaced face counts are invalid or truncated in %1, part %2.")
                                            .arg(geometryPath)
                                            .arg(partId);
                    }
                    return false;
                }

                qint64 connectivityCount = 0;
                qint32 minFaceNodes = 0;
                qint32 maxFaceNodes = 0;
                if (!sumEnSightGoldInt32Block(&file, faceCount, &connectivityCount, &minFaceNodes, &maxFaceNodes)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QString("EnSight Gold nfaced node counts are invalid or truncated in %1, part %2.")
                                            .arg(geometryPath)
                                            .arg(partId);
                    }
                    return false;
                }

                if (faceCount > std::numeric_limits<qint32>::max()
                    || connectivityCount > std::numeric_limits<qint32>::max()) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QString("EnSight Gold geometry is valid but too large for the current VTK EnSight reader. "
                                                "Part %1 contains %2 nfaced cells, %3 faces, and %4 face-node connectivity indices. "
                                                "The connectivity index count exceeds the 32-bit limit (%5), which triggers VTK's "
                                                "\"bad array new length\" failure during UpdateTimeStep. Re-export the STAR-CCM+ data "
                                                "as partitioned/split EnSight parts, VTU/VTK, CGNS, or another chunked format before importing.")
                                            .arg(partId)
                                            .arg(elementCount)
                                            .arg(faceCount)
                                            .arg(connectivityCount)
                                            .arg(std::numeric_limits<qint32>::max());
                    }
                    return false;
                }

                if (!skipEnSightGoldBytes(&file, connectivityCount * 4)) {
                    if (errorMessage != nullptr) {
                        *errorMessage = QString("EnSight Gold nfaced connectivity is truncated in %1, part %2.")
                                            .arg(geometryPath)
                                            .arg(partId);
                    }
                    return false;
                }
                continue;
            }

            const int nodesPerElement = enSightGoldFixedNodesPerElement(text);
            if (nodesPerElement <= 0) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("Unsupported EnSight Gold element type \"%1\" in %2, part %3.")
                                        .arg(text, geometryPath)
                                        .arg(partId);
                }
                return false;
            }
            if (!skipEnSightGoldBytes(&file, static_cast<qint64>(elementCount) * nodesPerElement * 4)) {
                if (errorMessage != nullptr) {
                    *errorMessage = QString("EnSight Gold connectivity is truncated in %1, part %2, element type %3.")
                                        .arg(geometryPath)
                                        .arg(partId)
                                        .arg(text);
                }
                return false;
            }
        }
    }
    return true;
}

bool checkEnSightGoldSourcesForVtk(const QVector<EnSightGoldCaseSource>& cases, QString* errorMessage) {
    QStringList checkedGeometryFiles;
    for (const EnSightGoldCaseSource& source : cases) {
        for (const QString& pattern : source.geometryPatterns) {
            const QString geometryPath = resolveEnSightGoldReferencedPath(source, pattern, 0);
            if (checkedGeometryFiles.contains(geometryPath, Qt::CaseInsensitive)) {
                continue;
            }
            checkedGeometryFiles.append(geometryPath);
            if (!checkEnSightGoldGeometryForVtk(geometryPath, errorMessage)) {
                return false;
            }
        }
    }
    return true;
}

bool timeCodesEquivalent(double lhs, double rhs) {
    return std::abs(lhs - rhs) <= 1e-12;
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
        cases->append(source);
    }

    for (int caseIndex = 0; caseIndex < cases->size(); ++caseIndex) {
        const EnSightGoldCaseSource& source = cases->at(caseIndex);
        for (const double timeValue : source.timeValues) {
            frames->append(EnSightGoldFrameSource{caseIndex, timeValue});
        }
    }

    std::sort(frames->begin(), frames->end(), [](const EnSightGoldFrameSource& lhs, const EnSightGoldFrameSource& rhs) {
        if (!timeCodesEquivalent(lhs.timeValue, rhs.timeValue)) {
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

}  // namespace

FlowDataScanResult FlowDataInspector::scanFolder(const QString& folder,
                                                 FlowDataKind dataKind,
                                                 bool coordinateCheck,
                                                 const ProgressCallback& progress,
                                                 const CancellationCallback& shouldCancel) {
    FlowDataScanResult result;
    result.sourceFolder = normalizedPath(folder);
    result.entry.dataKind = dataKind;
    result.entry.source = FlowDataSource::Link;
    result.entry.dataName = QFileInfo(folder).fileName();
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

    const QDir dir(folder);
    if (!dir.exists()) {
        result.errorMessage = QString("Flow data folder does not exist: %1").arg(folder);
        return result;
    }

    result.files = collectVtkFiles(folder, dataKind == FlowDataKind::Raw);
    if (result.files.isEmpty()) {
        const QString extensions = dataKind == FlowDataKind::Raw
            ? QStringLiteral("*.vtk, *.vti, *.vtr, *.vts, *.vtu, *.vtm, *.vtkhdf, *.nek5000, *.nek")
            : QStringLiteral("*.vtk, *.vti, *.vtr, *.vts, *.vtu, *.vtm, *.vtkhdf");
        result.errorMessage = QString("No VTK series files (%1) were found in %2.").arg(extensions, folder);
        return result;
    }
    if (std::all_of(result.files.cbegin(), result.files.cend(), [](const QString& file) {
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
    if (progress) {
        progress(QStringLiteral("Validating flow data"), 1, coordinateCheck ? result.files.size() : 2);
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
    if (dataKind == FlowDataKind::Resampled && !first.uniformStructured) {
        result.errorMessage = "Resampled flow data must be a uniform structured VTK dataset.";
        return result;
    }

    const double* bounds = first.bounds;
    result.entry.xRange = formatRange(bounds[0], bounds[1]);
    result.entry.yRange = formatRange(bounds[2], bounds[3]);
    result.entry.zRange = formatRange(bounds[4], bounds[5]);
    const bool reducedZ = (dims[2] > 0 && dims[2] <= 1) || std::abs(bounds[5] - bounds[4]) <= 1e-12;
    result.entry.dimension = reducedZ ? QStringLiteral("2") : QStringLiteral("3");
    const int dimension = result.entry.dimension == QStringLiteral("2") ? 2 : 3;

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
        result.errorMessage = "The first VTK snapshot does not contain numeric point or cell arrays.";
        return result;
    }

    if (coordinateCheck) {
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
            if (!boundsEquivalent(first, frame)) {
                result.errorMessage = QString("VTK snapshot grid coordinates differ from the first frame: %1")
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

FlowDataScanResult FlowDataInspector::scanNek5000File(const QString& filePath,
                                                      FlowDataKind dataKind,
                                                      bool /*coordinateCheck*/,
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

    if (progress) {
        progress(QStringLiteral("Validating Nek5000 flow data"), 1, 2);
    }
    if (canceled()) {
        return result;
    }

    QString readError;
    const QVector<double> timeValues = nek5000TimeValues(normalizedFilePath, &readError);
    if (!readError.trimmed().isEmpty() && timeValues.isEmpty()) {
        result.errorMessage = readError;
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

    vtkSmartPointer<vtkDataSet> first = readAnyVtk(normalizedFilePath, &readError);
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
    const bool reducedZ = (dims[2] > 0 && dims[2] <= 1) || std::abs(bounds[5] - bounds[4]) <= 1e-12;
    result.entry.dimension = reducedZ ? QStringLiteral("2") : QStringLiteral("3");
    const int dimension = result.entry.dimension == QStringLiteral("2") ? 2 : 3;

    if (dataKind == FlowDataKind::Resampled) {
        result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                     .arg(std::max(1, dims[0]))
                                     .arg(std::max(1, dims[1]))
                                     .arg(std::max(1, dims[2]));
    } else {
        result.entry.meshCount = QString::number(static_cast<qlonglong>(std::max<qsizetype>(1, timeValues.size())));
    }

    QList<ArrayInfo> arrays;
    collectArrays(first->GetPointData(), &arrays, &result.termCandidates);
    collectArrays(first->GetCellData(), &arrays, &result.termCandidates);
    result.termCandidates.sort(Qt::CaseInsensitive);
    inferTerms(arrays, dimension, &result.entry);
    result.entry.termCandidates = result.termCandidates;

    if (result.termCandidates.isEmpty()) {
        result.errorMessage = "The Nek5000 snapshot does not contain numeric point or cell arrays.";
        return result;
    }

    if (progress) {
        progress(QStringLiteral("Validating Nek5000 flow data"), 2, 2);
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
    const bool reducedZ = (dims[2] > 0 && dims[2] <= 1) || std::abs(bounds[5] - bounds[4]) <= 1e-12;
    result.entry.dimension = reducedZ ? QStringLiteral("2") : QStringLiteral("3");
    const int dimension = result.entry.dimension == QStringLiteral("2") ? 2 : 3;

    if (dataKind == FlowDataKind::Resampled) {
        result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                     .arg(std::max(1, dims[0]))
                                     .arg(std::max(1, dims[1]))
                                     .arg(std::max(1, dims[2]));
    } else {
        result.entry.meshCount = QString::number(timeValues.size());
    }

    QList<ArrayInfo> arrays;
    collectArrays(first->GetPointData(), &arrays, &result.termCandidates);
    collectArrays(first->GetCellData(), &arrays, &result.termCandidates);
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
            if (coordinateCheck && !boundsEquivalent(first, frame)) {
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
    result.sourceFolder = casePaths.size() == 1 ? normalizedPath(casePaths.first()) : casePaths.join(QStringLiteral("; "));
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
    result.files.reserve(cases.size());
    for (const EnSightGoldCaseSource& source : std::as_const(cases)) {
        result.files.append(source.path);
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
    result.entry.dt = formatDouble(inferredTimeInterval(timeValues));

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
        inferTerms(arrays, 3, &result.entry);
        result.entry.termCandidates = result.termCandidates;
        result.entry.meshCount = QString::number(frames.size());
        result.entry.dimension = QStringLiteral("3");
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
        firstBoundsText = boundsSummary(firstBounds);
        result.entry.xRange = formatRange(firstBounds[0], firstBounds[1]);
        result.entry.yRange = formatRange(firstBounds[2], firstBounds[3]);
        result.entry.zRange = formatRange(firstBounds[4], firstBounds[5]);
        const bool reducedZ = (dims[2] > 0 && dims[2] <= 1) || std::abs(firstBounds[5] - firstBounds[4]) <= 1e-12;
        result.entry.dimension = reducedZ ? QStringLiteral("2") : QStringLiteral("3");
        dimension = result.entry.dimension == QStringLiteral("2") ? 2 : 3;

        if (dataKind == FlowDataKind::Resampled) {
            result.entry.meshCount = QStringLiteral("%1 %2 %3")
                                         .arg(std::max(1, dims[0]))
                                         .arg(std::max(1, dims[1]))
                                         .arg(std::max(1, dims[2]));
        } else {
            result.entry.meshCount = QString::number(frames.size());
        }

        collectArrays(first->GetPointData(), &arrays, &result.termCandidates);
        collectArrays(first->GetCellData(), &arrays, &result.termCandidates);
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
            if (coordinateCheck && !boundsRangeEquivalent(firstBounds, frame) && result.warningMessage.isEmpty()) {
                result.warningMessage = QString("EnSight Gold frame XYZ bounds differ from the first frame at time %1 in %2.\nFirst bounds: %3\nCurrent bounds: %4\nImport will continue; verify the combined dataset before analysis.")
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
