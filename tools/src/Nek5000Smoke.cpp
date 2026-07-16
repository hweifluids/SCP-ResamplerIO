#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkErrorCode.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkNek5000Reader.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkResampleToImage.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>

#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QSaveFile>
#include <QString>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>

#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr std::array<int, 3> kDimensions{250, 80, 50};
constexpr std::array<double, 6> kBounds{
    -3.470588207244873,
    12.45714282989502,
    -2.858823537826538,
    2.858823537826538,
    0.0,
    3.418181896209717};
constexpr int kFirstFrame = 490;
constexpr int kLastFrame = 509;
constexpr int kFrameCount = kLastFrame - kFirstFrame + 1;
constexpr vtkIdType kPointCount =
    static_cast<vtkIdType>(kDimensions[0]) * kDimensions[1] * kDimensions[2];
constexpr double kGeometryTolerance = 1.0e-6;
constexpr double kTimeTolerance = 1.0e-6;

QString toQString(const std::filesystem::path& path)
{
  return QString::fromUtf8(path.u8string().c_str());
}

std::string frameFileName(int frame)
{
  std::ostringstream name;
  name << "hemisphere_re1370_" << std::setw(4) << std::setfill('0') << frame << ".vti";
  return name.str();
}

void require(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

bool nearlyEqual(double actual, double expected, double tolerance)
{
  return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
}

void validateBounds(const double actual[6], const std::string& context)
{
  for (std::size_t axis = 0; axis < kBounds.size(); ++axis)
  {
    if (!nearlyEqual(actual[axis], kBounds[axis], kGeometryTolerance))
    {
      std::ostringstream message;
      message << context << " bound[" << axis << "] is " << std::setprecision(17)
              << actual[axis] << ", expected " << kBounds[axis];
      throw std::runtime_error(message.str());
    }
  }
}

std::vector<double> readAvailableTimes(vtkNek5000Reader* reader)
{
  reader->UpdateInformation();
  vtkInformation* information = reader->GetOutputInformation(0);
  auto* timeStepsKey = vtkStreamingDemandDrivenPipeline::TIME_STEPS();
  require(information != nullptr && information->Has(timeStepsKey),
          "The Nek5000 descriptor does not advertise time steps.");

  const int count = information->Length(timeStepsKey);
  std::vector<double> times;
  times.reserve(static_cast<std::size_t>(count));
  for (int index = 0; index < count; ++index)
  {
    times.push_back(information->Get(timeStepsKey, index));
  }
  return times;
}

void addTimeValue(vtkImageData* image, double time)
{
  image->GetFieldData()->RemoveArray("TimeValue");
  vtkNew<vtkDoubleArray> timeValue;
  timeValue->SetName("TimeValue");
  timeValue->SetNumberOfTuples(1);
  timeValue->SetValue(0, time);
  image->GetFieldData()->AddArray(timeValue);
}

void writePvd(const std::filesystem::path& outputDirectory, const std::vector<double>& times)
{
  const auto path = outputDirectory / "hemisphere_re1370_250x80x50.pvd";
  QSaveFile file(toQString(path));
  require(file.open(QIODevice::WriteOnly | QIODevice::Text),
          "Cannot open the PVD collection for writing: " + path.string());

  QXmlStreamWriter xml(&file);
  xml.setAutoFormatting(true);
  xml.writeStartDocument();
  xml.writeStartElement(QStringLiteral("VTKFile"));
  xml.writeAttribute(QStringLiteral("type"), QStringLiteral("Collection"));
  xml.writeAttribute(QStringLiteral("version"), QStringLiteral("0.1"));
  xml.writeAttribute(QStringLiteral("byte_order"), QStringLiteral("LittleEndian"));
  xml.writeStartElement(QStringLiteral("Collection"));
  for (int offset = 0; offset < kFrameCount; ++offset)
  {
    xml.writeEmptyElement(QStringLiteral("DataSet"));
    xml.writeAttribute(QStringLiteral("timestep"),
                       QString::number(times[static_cast<std::size_t>(offset)], 'g', 17));
    xml.writeAttribute(QStringLiteral("group"), QString());
    xml.writeAttribute(QStringLiteral("part"), QStringLiteral("0"));
    xml.writeAttribute(QStringLiteral("file"),
                       QString::fromStdString(frameFileName(kFirstFrame + offset)));
  }
  xml.writeEndElement();
  xml.writeEndElement();
  xml.writeEndDocument();
  require(!xml.hasError() && file.commit(), "Failed to commit the PVD collection: " + path.string());
}

void writeProvenance(const std::filesystem::path& descriptor,
                     const std::filesystem::path& outputDirectory,
                     const std::vector<double>& times)
{
  QJsonArray dimensions;
  for (int value : kDimensions)
  {
    dimensions.append(value);
  }
  QJsonArray bounds;
  for (double value : kBounds)
  {
    bounds.append(value);
  }
  QJsonArray frames;
  for (int offset = 0; offset < kFrameCount; ++offset)
  {
    QJsonObject frame;
    frame.insert(QStringLiteral("source_index"), kFirstFrame + offset);
    frame.insert(QStringLiteral("time"), times[static_cast<std::size_t>(offset)]);
    frame.insert(QStringLiteral("file"),
                 QString::fromStdString(frameFileName(kFirstFrame + offset)));
    frames.append(frame);
  }

  QJsonObject velocity;
  velocity.insert(QStringLiteral("name"), QStringLiteral("Velocity"));
  velocity.insert(QStringLiteral("association"), QStringLiteral("point"));
  velocity.insert(QStringLiteral("components"), 3);
  velocity.insert(QStringLiteral("scalar_type"), QStringLiteral("float32"));

  QJsonObject root;
  root.insert(QStringLiteral("schema"), QStringLiteral("scp-resamplerio-smoke-v1"));
  root.insert(QStringLiteral("source_descriptor"),
              QString::fromStdString(descriptor.filename().string()));
  root.insert(QStringLiteral("selection"), QStringLiteral("middle-20-of-1000"));
  root.insert(QStringLiteral("dimensions"), dimensions);
  root.insert(QStringLiteral("bounds"), bounds);
  root.insert(QStringLiteral("field"), velocity);
  root.insert(QStringLiteral("frames"), frames);

  const auto path = outputDirectory / "dataset.json";
  QSaveFile file(toQString(path));
  require(file.open(QIODevice::WriteOnly | QIODevice::Text),
          "Cannot open provenance JSON for writing: " + path.string());
  require(file.write(QJsonDocument(root).toJson(QJsonDocument::Indented)) >= 0 && file.commit(),
          "Failed to commit provenance JSON: " + path.string());
}

void generate(const std::filesystem::path& descriptor,
              const std::filesystem::path& outputDirectory)
{
  require(std::filesystem::is_regular_file(descriptor),
          "Nek5000 descriptor does not exist: " + descriptor.string());
  require(descriptor.extension() == ".nek5000",
          "The generator input must be a .nek5000 descriptor.");
  std::filesystem::create_directories(outputDirectory);

  vtkNew<vtkNek5000Reader> reader;
  reader->SetFileName(descriptor.string().c_str());
  reader->SetCleanGrid(0);
  reader->EnableAllPointArrays();
  const std::vector<double> allTimes = readAvailableTimes(reader);
  require(allTimes.size() > static_cast<std::size_t>(kLastFrame),
          "The Nek5000 sequence does not contain source frames 490 through 509.");

  std::vector<double> selectedTimes;
  selectedTimes.reserve(kFrameCount);
  for (int frame = kFirstFrame; frame <= kLastFrame; ++frame)
  {
    const double time = allTimes[static_cast<std::size_t>(frame)];
    std::cout << "Resampling source frame " << frame << " at time "
              << std::setprecision(17) << time << std::endl;
    require(reader->UpdateTimeStep(time) != 0,
            "VTK failed to update Nek5000 source frame " + std::to_string(frame) + ".");

    vtkUnstructuredGrid* source = reader->GetOutput();
    require(source != nullptr && source->GetNumberOfPoints() > 0,
            "Nek5000 source frame is empty: " + std::to_string(frame));
    vtkDataArray* sourceVelocity = source->GetPointData()->GetArray("Velocity");
    require(sourceVelocity != nullptr && sourceVelocity->GetNumberOfComponents() == 3,
            "Nek5000 source frame does not contain point-associated Velocity with 3 components.");
    double sourceBounds[6];
    source->GetBounds(sourceBounds);
    validateBounds(sourceBounds, "Nek5000 source frame " + std::to_string(frame));

    vtkNew<vtkResampleToImage> resampler;
    resampler->SetInputDataObject(source);
    resampler->SetUseInputBounds(false);
    resampler->SetSamplingDimensions(kDimensions[0], kDimensions[1], kDimensions[2]);
    resampler->SetSamplingBounds(kBounds[0], kBounds[1], kBounds[2], kBounds[3],
                                 kBounds[4], kBounds[5]);
    resampler->Update();

    vtkImageData* image = resampler->GetOutput();
    require(image != nullptr && image->GetNumberOfPoints() == kPointCount,
            "Resampling did not produce the required 250 x 80 x 50 image.");
    addTimeValue(image, time);

    const auto outputPath = outputDirectory / frameFileName(frame);
    vtkNew<vtkXMLImageDataWriter> writer;
    writer->SetFileName(outputPath.string().c_str());
    writer->SetInputData(image);
    writer->SetDataModeToAppended();
    writer->EncodeAppendedDataOn();
    writer->SetCompressorTypeToZLib();
    require(writer->Write() == 1 && writer->GetErrorCode() == vtkErrorCode::NoError,
            "Failed to write sampled frame: " + outputPath.string());
    selectedTimes.push_back(time);
  }

  writePvd(outputDirectory, selectedTimes);
  writeProvenance(descriptor, outputDirectory, selectedTimes);
  std::cout << "Generated " << kFrameCount << " frames in " << outputDirectory.string()
            << std::endl;
}

double validateFrame(const std::filesystem::path& path, int frame)
{
  require(std::filesystem::is_regular_file(path), "Missing smoke frame: " + path.string());
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(path.string().c_str());
  reader->Update();
  require(reader->GetErrorCode() == vtkErrorCode::NoError,
          "VTK failed to read smoke frame: " + path.string());

  vtkImageData* image = reader->GetOutput();
  require(image != nullptr, "Smoke frame has no vtkImageData output: " + path.string());
  int dimensions[3];
  image->GetDimensions(dimensions);
  for (std::size_t axis = 0; axis < kDimensions.size(); ++axis)
  {
    require(dimensions[axis] == kDimensions[axis],
            "Smoke frame " + std::to_string(frame) + " has incorrect dimensions.");
  }
  require(image->GetNumberOfPoints() == kPointCount,
          "Smoke frame " + std::to_string(frame) + " has an incorrect point count.");
  double bounds[6];
  image->GetBounds(bounds);
  validateBounds(bounds, "Smoke frame " + std::to_string(frame));

  vtkDataArray* velocity = image->GetPointData()->GetArray("Velocity");
  require(velocity != nullptr && velocity->GetNumberOfComponents() == 3 &&
              velocity->GetNumberOfTuples() == kPointCount && velocity->GetDataType() == VTK_FLOAT,
          "Smoke frame " + std::to_string(frame) +
              " must contain point-associated Velocity(float3).");

  vtkDataArray* validMask = image->GetPointData()->GetArray("vtkValidPointMask");
  require(validMask != nullptr && validMask->GetNumberOfTuples() == kPointCount,
          "Smoke frame " + std::to_string(frame) + " is missing vtkValidPointMask.");
  double maskRange[2];
  validMask->GetRange(maskRange);
  require(maskRange[0] >= 0.0 && maskRange[1] <= 1.0,
          "Smoke frame " + std::to_string(frame) + " has a non-binary valid-point mask.");
  vtkIdType validPointCount = 0;
  for (vtkIdType point = 0; point < kPointCount; ++point)
  {
    validPointCount += validMask->GetTuple1(point) > 0.5 ? 1 : 0;
  }
  require(validPointCount * 100 >= kPointCount * 99,
          "Smoke frame " + std::to_string(frame) +
              " has less than 99 percent valid samples over the full source bounds.");
  std::cout << "  valid samples: " << validPointCount << "/" << kPointCount << std::endl;

  vtkDataArray* timeValue = image->GetFieldData()->GetArray("TimeValue");
  require(timeValue != nullptr && timeValue->GetNumberOfTuples() == 1,
          "Smoke frame " + std::to_string(frame) + " is missing TimeValue field data.");
  return timeValue->GetComponent(0, 0);
}

void validatePvd(const std::filesystem::path& outputDirectory,
                 const std::vector<double>& frameTimes)
{
  const auto path = outputDirectory / "hemisphere_re1370_250x80x50.pvd";
  QFile file(toQString(path));
  require(file.open(QIODevice::ReadOnly | QIODevice::Text),
          "Missing or unreadable PVD collection: " + path.string());
  QXmlStreamReader xml(&file);
  int dataSetIndex = 0;
  while (!xml.atEnd())
  {
    xml.readNext();
    if (!xml.isStartElement() || xml.name() != QStringLiteral("DataSet"))
    {
      continue;
    }
    require(dataSetIndex < kFrameCount, "PVD contains more than 20 DataSet entries.");
    const auto attributes = xml.attributes();
    const QString expectedFile =
        QString::fromStdString(frameFileName(kFirstFrame + dataSetIndex));
    require(attributes.value(QStringLiteral("file")) == expectedFile,
            "PVD frame ordering or filename is incorrect.");
    bool timeOk = false;
    const double time = attributes.value(QStringLiteral("timestep")).toDouble(&timeOk);
    require(timeOk && nearlyEqual(time, frameTimes[static_cast<std::size_t>(dataSetIndex)],
                                  kTimeTolerance),
            "PVD timestep does not match the VTI TimeValue.");
    ++dataSetIndex;
  }
  require(!xml.hasError(), "PVD XML is malformed: " + xml.errorString().toStdString());
  require(dataSetIndex == kFrameCount, "PVD must contain exactly 20 DataSet entries.");
}

void validateProvenance(const std::filesystem::path& outputDirectory)
{
  const auto path = outputDirectory / "dataset.json";
  QFile file(toQString(path));
  require(file.open(QIODevice::ReadOnly | QIODevice::Text),
          "Missing or unreadable provenance JSON: " + path.string());
  QJsonParseError error;
  const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &error);
  require(error.error == QJsonParseError::NoError && document.isObject() &&
              !document.object().isEmpty(),
          "dataset.json must be a non-empty JSON object: " + error.errorString().toStdString());
}

void validate(const std::filesystem::path& outputDirectory)
{
  require(std::filesystem::is_directory(outputDirectory),
          "Smoke output directory does not exist: " + outputDirectory.string());
  std::vector<double> frameTimes;
  frameTimes.reserve(kFrameCount);
  for (int frame = kFirstFrame; frame <= kLastFrame; ++frame)
  {
    const auto path = outputDirectory / frameFileName(frame);
    const double time = validateFrame(path, frame);
    if (!frameTimes.empty())
    {
      require(time > frameTimes.back(), "Smoke frame TimeValue values must be strictly increasing.");
    }
    frameTimes.push_back(time);
    std::cout << "Validated " << path.filename().string() << std::endl;
  }
  require(nearlyEqual(frameTimes.front(), 524.5499999911, kTimeTolerance),
          "First smoke time must correspond to source frame 490.");
  require(nearlyEqual(frameTimes.back(), 525.499999991, kTimeTolerance),
          "Last smoke time must correspond to source frame 509.");
  validatePvd(outputDirectory, frameTimes);
  validateProvenance(outputDirectory);
  std::cout << "Smoke dataset validation passed: 20 frames, 250 x 80 x 50, Velocity(float3)."
            << std::endl;
}

void printUsage(const char* program)
{
  std::cout << "Usage:\n"
            << "  " << program << " generate <hemi.nek5000> <output-directory>\n"
            << "  " << program << " validate <output-directory>\n";
}

} // namespace

int main(int argc, char* argv[])
{
  if (argc == 2 && std::string(argv[1]) == "--solver-info")
  {
    return 2;
  }
  if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))
  {
    printUsage(argv[0]);
    return 0;
  }

  try
  {
    if (argc == 4 && std::string(argv[1]) == "generate")
    {
      generate(argv[2], argv[3]);
      return 0;
    }
    if (argc == 3 && std::string(argv[1]) == "validate")
    {
      validate(argv[2]);
      return 0;
    }
    printUsage(argv[0]);
    return 2;
  }
  catch (const std::exception& error)
  {
    std::cerr << "SCP Nek5000 smoke error: " << error.what() << std::endl;
    return 1;
  }
}
