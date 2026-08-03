#pragma once

#include <QString>
#include <QStringList>
#include <QVector>

enum class FlowDataKind {
    Raw,
    Resampled
};

enum class FlowDataSource {
    Link,
    Import
};

struct FlowDataSourceFingerprint {
    QString path;
    qint64 size = 0;
    QString sha256Hex;

    bool operator==(const FlowDataSourceFingerprint& other) const {
        return path == other.path && size == other.size && sha256Hex == other.sha256Hex;
    }
};

struct FlowDataEntry {
    QString id;
    FlowDataKind dataKind = FlowDataKind::Raw;
    FlowDataSource source = FlowDataSource::Link;
    QString dataName;
    QString folder;
    QString dt;
    QString xRange;
    QString yRange;
    QString zRange;
    QString uTerm;
    QString vTerm;
    QString wTerm;
    QString pTerm;
    QString dataVersion;
    QString meshCount;
    QString dimension;
    QString fileType = "vtk series";
    QString checkStatus = "unchecked";
    QString scpfPath;
    QStringList termCandidates;
};
