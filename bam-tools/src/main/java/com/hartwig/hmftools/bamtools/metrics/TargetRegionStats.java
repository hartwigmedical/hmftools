package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.DUAL_STRAND_COLUMN;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.DUPLICATES_COLUMN;
import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.TOTAL_READS_COLUMN;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.metrics.ImmutableTargetRegionMetrics;
import com.hartwig.hmftools.common.metrics.TargetRegionMetrics;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetRegionStats
{
    public final ChrBaseRegion Region;
    public ReadCounts Counts;

    public TargetRegionStats(final ChrBaseRegion region)
    {
        Region = region;
        Counts = new ReadCounts();
    }

    public String toString()
    {
        return format("%s reads(%s)", Region, Counts);
    }

    public TargetRegionMetrics convert()
    {
        return ImmutableTargetRegionMetrics.builder()
                .region(Region)
                .totalReads(Counts.TotalReads)
                .duplicateReads(Counts.Duplicates)
                .dualStrandReads(Counts.DualStrand)
                .build();
    }
}
