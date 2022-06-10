package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class DisruptionData implements ComparableItem
{
    public final StructuralVariantData SvData;
    public final LinxBreakend Breakend;

    protected static final String FLD_UNDISRUPTED_CN = "undisruptedCopyNumber";

    public DisruptionData(final StructuralVariantData svData, final LinxBreakend breakend)
    {
        SvData = svData;
        Breakend = breakend;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public String key()
    {
        return String.format("%s %d_%s %s_%d - %s_%d %s",
                SvData.id(), SvData.type(), SvData.startChromosome(), SvData.startPosition(), SvData.endChromosome(), SvData.endPosition(),
                Breakend.gene());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
//        values.add(String.format("Qual(%.0f)", Variant.qual()));
//        values.add(String.format("Tier(%s)", Variant.tier().toString()));
//        values.add(String.format("TotalReadCount(%d)", Variant.totalReadCount()));
//        values.add(String.format("AlleleReadCount(%d)", Variant.alleleReadCount()));
        return values;
    }

    @Override
    public boolean reportable() { return Breakend.reportedDisruption(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherSv = (DisruptionData)other;

        if(otherSv.SvData.type() != SvData.type())
            return false;

        if(!otherSv.SvData.startChromosome().equals(SvData.startChromosome()) || !otherSv.SvData.endChromosome().equals(SvData.endChromosome()))
            return false;

        if(otherSv.SvData.startPosition() != SvData.startPosition() || otherSv.SvData.endPosition() != SvData.endPosition())
            return false;

        if(otherSv.SvData.startOrientation() != SvData.startOrientation() || otherSv.SvData.endOrientation() != SvData.endOrientation())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", reportable(), otherBreakend.reportable());
        checkDiff(diffs, "geneOrientation", Breakend.geneOrientation(), otherBreakend.Breakend.geneOrientation());
        checkDiff(diffs, "nextSpliceExonRank", Breakend.nextSpliceExonRank(), otherBreakend.Breakend.nextSpliceExonRank());
        checkDiff(diffs, "regionType", Breakend.regionType(), otherBreakend.Breakend.regionType());
        checkDiff(diffs, "codingContext", Breakend.codingContext(), otherBreakend.Breakend.codingContext());

        checkDiff(diffs, FLD_UNDISRUPTED_CN, Breakend.undisruptedCopyNumber(), otherBreakend.Breakend.undisruptedCopyNumber(), thresholds);

        if(matchLevel == REPORTABLE)
            return null;


        return null;
    }
}
