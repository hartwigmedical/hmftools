package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

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

    protected static final String FLD_REGION_TYPE = "RegionType";
    protected static final String FLD_CODING_CONTEXT = "CodingContext";
    protected static final String FLD_GENE_ORIENT = "GeneOrientation";
    protected static final String FLD_NEXT_SPLICE = "NextSpliceExonRank";

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
        return String.format("%s %d_%s %s:%d-%s:%d",
                Breakend.gene(), SvData.id(), SvData.type(),
                SvData.startChromosome(), SvData.startPosition(), SvData.endChromosome(), SvData.endPosition());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("%s", Breakend.reportedDisruption()));
        values.add(String.format("%s", Breakend.regionType()));
        values.add(String.format("%s", Breakend.codingContext()));
        values.add(String.format("%s", Breakend.geneOrientation()));
        values.add(String.format("%d", Breakend.nextSpliceExonRank()));
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

        checkDiff(diffs, FLD_REGION_TYPE, Breakend.regionType(), otherBreakend.Breakend.regionType());
        checkDiff(diffs, FLD_CODING_CONTEXT, Breakend.codingContext(), otherBreakend.Breakend.codingContext());
        checkDiff(diffs, FLD_REPORTED, reportable(), otherBreakend.reportable());
        checkDiff(diffs, FLD_GENE_ORIENT, Breakend.geneOrientation(), otherBreakend.Breakend.geneOrientation());
        checkDiff(diffs, FLD_NEXT_SPLICE, Breakend.nextSpliceExonRank(), otherBreakend.Breakend.nextSpliceExonRank());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
