package com.hartwig.hmftools.common.purple;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

public class GeneCopyNumber implements GenomeRegion
{
    public final String Chromosome;
    public final int PositionStart;
    public final int PositionEnd;
    public final String GeneName;
    public final String TransName;
    public final boolean IsCanonical;

    public final String ChromosomeBand;

    public final double MaxCopyNumber;
    public final double MinCopyNumber;
    public final double MinMinorAlleleCopyNumber;
    public final double RelativeMinCopyNumber;

    public final int SomaticRegions;
    public final int MinRegions;
    public final int MinRegionStart;
    public final int MinRegionEnd;

    public final int DepthWindowCount;
    public final double GcContent;

    public final SegmentSupport MinRegionStartSupport;
    public final SegmentSupport MinRegionEndSupport;
    public final CopyNumberMethod MinRegionMethod;

    private DriverType mDriverType;
    private ReportedStatus mReportedStatus;

    public GeneCopyNumber(
            final String chromosome, final int positionStart, final int positionEnd, final String geneName, final String transName,
            final boolean isCanonical, final String chromosomeBand, final double maxCopyNumber, final double minCopyNumber,
            final double minMinorAlleleCopyNumber, final int somaticRegions, final int minRegions, final int minRegionStart,
            final int minRegionEnd, final int depthWindowCount, final double gcContent, final SegmentSupport minRegionStartSupport,
            final SegmentSupport minRegionEndSupport, final CopyNumberMethod minRegionMethod, final double relativeMinCopyNumber)
    {
        Chromosome = chromosome;
        PositionStart = positionStart;
        PositionEnd = positionEnd;
        GeneName = geneName;
        TransName = transName;
        IsCanonical = isCanonical;
        ChromosomeBand = chromosomeBand;
        MaxCopyNumber = maxCopyNumber;
        MinCopyNumber = minCopyNumber;
        MinMinorAlleleCopyNumber = minMinorAlleleCopyNumber;
        SomaticRegions = somaticRegions;
        MinRegions = minRegions;
        MinRegionStart = minRegionStart;
        MinRegionEnd = minRegionEnd;
        DepthWindowCount = depthWindowCount;
        GcContent = gcContent;
        MinRegionStartSupport = minRegionStartSupport;
        MinRegionEndSupport = minRegionEndSupport;
        MinRegionMethod = minRegionMethod;
        RelativeMinCopyNumber = relativeMinCopyNumber;

        mDriverType = DriverType.UNKNOWN;
        mReportedStatus = ReportedStatus.NONE;
    }

    public String chromosome() { return Chromosome; }
    public int start() { return PositionStart; }
    public int end() { return PositionEnd; }

    // convenience for existing code
    public String geneName() { return GeneName; }
    public double maxCopyNumber() { return MaxCopyNumber; }
    public double minCopyNumber() { return MinCopyNumber; }

    public int minRegionBases()
    {
        return MinRegionEnd - MinRegionStart + 1;
    }

    public void setDriverType(final DriverType type) { mDriverType = type; }
    public DriverType driverType() { return mDriverType; }

    public ReportedStatus reportedStatus() { return mReportedStatus; }
    public void setReportedStatus(final ReportedStatus status) { mReportedStatus = status; }

    public int totalRegions()
    {
        return SomaticRegions;
    }
}
