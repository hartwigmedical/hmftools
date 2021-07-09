package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

import org.apache.commons.compress.utils.Lists;

public class DisruptionData implements ComparableItem
{
    public final LinxBreakend Breakend;

    public DisruptionData(final LinxBreakend breakend)
    {
        Breakend = breakend;
    }

    public Category category() { return DISRUPTION; }

    public boolean reportable() { return Breakend.reportedDisruption(); }

    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        if(!otherBreakend.Breakend.gene().equals(Breakend.gene()))
            return false;

        if(!otherBreakend.Breakend.transcriptId().equals(Breakend.transcriptId()))
            return false;

        if(otherBreakend.Breakend.nextSpliceExonRank() != Breakend.nextSpliceExonRank())
            return false;

        if(!otherBreakend.Breakend.codingContext().equals(Breakend.codingContext()))
            return false;

        if(!otherBreakend.Breakend.regionType().equals(Breakend.regionType()))
            return false;

        if(otherBreakend.Breakend.isStart() != Breakend.isStart())
            return false;

        return true;
    }

    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Breakend.reportedDisruption(), otherBreakend.Breakend.reportedDisruption());
        checkDiff(diffs, "disruptive", Breakend.disruptive(), otherBreakend.Breakend.disruptive());

        if(matchLevel == REPORTABLE)
            return diffs;

        return diffs;
    }

    public String description()
    {
        return String.format("%d_%s_%s", Breakend.svId(), Breakend.isStart(), Breakend.gene());
    }
}
