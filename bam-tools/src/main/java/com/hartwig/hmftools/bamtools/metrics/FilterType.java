package com.hartwig.hmftools.bamtools.metrics;

import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_BASEQ_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_CAPPED_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_DUPE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_MAPQ_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_OVERLAP_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_UNPAIRED_COLUMN;

public enum FilterType
{
    UNFILTERED("Unfiltered", "UNFILTERED"),
    LOW_MAP_QUAL("LowMapQual", PCT_EXC_MAPQ_COLUMN),
    DUPLICATE("Duplicate", PCT_EXC_DUPE_COLUMN),
    MATE_UNMAPPED("MateUnmapped", PCT_EXC_UNPAIRED_COLUMN),
    LOW_BASE_QUAL("LowBaseQual", PCT_EXC_BASEQ_COLUMN),
    OVERLAPPED("Overlapped", PCT_EXC_OVERLAP_COLUMN),
    MAX_COVERAGE("MaxCoverage", PCT_EXC_CAPPED_COLUMN);

    private final String mDescription;
    private final String mFileColumn;

    FilterType(final String desc, final String fileColumn)
    {
        mDescription = desc;
        mFileColumn = fileColumn;
    }

    public String description() { return mDescription; }
    public String fileColumn() { return mFileColumn; }
}
