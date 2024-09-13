package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.compar.common.DiffThresholds;

public class LinxCommon
{
    protected static final String FLD_REGION_TYPE = "RegionType";
    protected static final String FLD_NEXT_SPLICE = "NextSpliceExonRank";
    protected static final String FLD_JUNCTION_COPY_NUMBER = "JunctionCopyNumber";
    protected static final String FLD_UNDISRUPTED_COPY_NUMBER = "UndisruptedCopyNumber";

    protected static List<String> comparedFieldNamesBreakends()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_REGION_TYPE, FLD_NEXT_SPLICE, FLD_JUNCTION_COPY_NUMBER, FLD_UNDISRUPTED_COPY_NUMBER);
    }

    protected static List<String> displayValuesBreakend(LinxBreakend breakend)
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("%s", breakend.reportedDisruption()));
        values.add(String.format("%s", breakend.regionType()));
        values.add(String.format("%d", breakend.nextSpliceExonRank()));
        values.add(String.format("%.2f", breakend.junctionCopyNumber()));
        values.add(String.format("%.2f", breakend.undisruptedCopyNumber()));
        return values;
    }

    protected static void checkDiffsBreakends(List<String> diffs, final LinxBreakend breakend, final LinxBreakend otherBreakend, final DiffThresholds thresholds)
    {
        checkDiff(diffs, FLD_REPORTED, breakend.reportedDisruption(), otherBreakend.reportedDisruption());
        checkDiff(diffs, FLD_REGION_TYPE, breakend.regionType().toString(), otherBreakend.regionType().toString());
        checkDiff(diffs, FLD_NEXT_SPLICE, breakend.nextSpliceExonRank(), otherBreakend.nextSpliceExonRank());
        checkDiff(diffs, FLD_JUNCTION_COPY_NUMBER, breakend.junctionCopyNumber(), otherBreakend.junctionCopyNumber(), thresholds);
        checkDiff(diffs, FLD_UNDISRUPTED_COPY_NUMBER, breakend.undisruptedCopyNumber(), otherBreakend.undisruptedCopyNumber(), thresholds);
        checkDiff(diffs, FLD_CHROMOSOME_BAND, breakend.chrBand(), otherBreakend.chrBand());
    }
}
