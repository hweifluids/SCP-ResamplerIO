#pragma once

#include <QtGlobal>

namespace Streamcenter::FlowDataProtocolLimits {

// A conversion request can contain both the user-selected case list and the
// materialized dependency inventory.  Keep their string budgets independent,
// then cap the complete JSON (including object/key overhead) as a second line
// of defence.
// Protocol strings reject C0 controls, so JSON escaping can expand the
// user-controlled portion by at most 2x (quotes/backslashes).  Sixteen MiB
// closes the worst case for the independent ordinary/inventory budgets plus
// all bounded fingerprint object, size, and key overhead without making the
// GUI synchronously parse an unnecessarily large response.
inline constexpr qint64 MaximumRequestBytes = 16 * 1024 * 1024;
inline constexpr qint64 MaximumResponseBytes = 16 * 1024 * 1024;
inline constexpr qint64 MaximumCaseMetadataBytes = 4 * 1024 * 1024;
inline constexpr qint64 MaximumCollectionMetadataBytes = 4 * 1024 * 1024;
inline constexpr qsizetype MaximumCaseCount = 4096;
inline constexpr qsizetype MaximumFieldCount = 4096;
inline constexpr qsizetype MaximumResultFileCount = 4096;
inline constexpr qsizetype MaximumTermCount = 4096;
inline constexpr qsizetype MaximumInventoryEntryCount = 16384;
inline constexpr qsizetype MaximumStringBytes = 64 * 1024;
inline constexpr qint64 MaximumAggregateStringBytes = 3 * 1024 * 1024;
inline constexpr qint64 MaximumResponseAggregateStringBytes = 7 * 1024 * 1024 / 2;
inline constexpr qint64 MaximumInventoryAggregateStringBytes = 4 * 1024 * 1024;

}  // namespace Streamcenter::FlowDataProtocolLimits
