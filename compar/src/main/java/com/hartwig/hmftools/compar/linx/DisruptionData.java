package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

import org.apache.commons.compress.utils.Lists;

public class DisruptionData implements ComparableItem
{
    public final String MappedGeneName;

    private final List<LinxBreakend> mBreakends;

    public DisruptionData(final String mappedName)
    {
        MappedGeneName = mappedName;
        mBreakends = Lists.newArrayList();
    }

    public final List<LinxBreakend> breakends() { return mBreakends; }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public boolean reportable() { return mBreakends.stream().anyMatch(x -> x.reportedDisruption()); }

    public boolean hasDisruptive() { return mBreakends.stream().anyMatch(x -> x.disruptive()); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        if(!otherBreakend.MappedGeneName.equals(MappedGeneName))
            return false;

        /*
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
        */

        return true;
    }

    @Override
    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", reportable(), otherBreakend.reportable());
        checkDiff(diffs, "disruptive", hasDisruptive(), otherBreakend.hasDisruptive());

        if(matchLevel == REPORTABLE)
            return diffs;

        return diffs;
    }

    @Override
    public String description()
    {
        return String.format("%s_%s_%d", MappedGeneName, reportable(), mBreakends.size());
        // return String.format("%d_%s_%s", Breakend.svId(), Breakend.isStart(), Breakend.gene());
    }

    @Override
    public String gene() { return MappedGeneName; }

}
