package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class DisruptionData implements ComparableItem
{
    public final StructuralVariantData SvData;
    public final LinxBreakend Breakend;
    private final GenomePositionImpl mComparisonStartGenomePosition;
    private final GenomePositionImpl mComparisonEndGenomePosition;
    private final boolean mCheckTranscript;

    protected static final String FLD_REGION_TYPE = "RegionType";
    protected static final String FLD_CODING_CONTEXT = "CodingContext";
    protected static final String FLD_GENE_ORIENT = "GeneOrientation";
    protected static final String FLD_NEXT_SPLICE = "NextSpliceExonRank";

    public DisruptionData(final StructuralVariantData svData, final LinxBreakend breakend,
            final GenomePositionImpl comparisonStartGenomePosition, final GenomePositionImpl comparisonEndGenomePosition,
            final boolean checkTranscript)
    {
        SvData = svData;
        Breakend = breakend;
        mComparisonStartGenomePosition = comparisonStartGenomePosition;
        mComparisonEndGenomePosition = comparisonEndGenomePosition;
        mCheckTranscript = checkTranscript;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public String key()
    {
        if(mComparisonStartGenomePosition.position() != SvData.startPosition() || mComparisonEndGenomePosition.position() != SvData.endPosition())
        {
            return String.format("%s %d_%s %s:%d-%s:%d liftover(%s:%d-%s:%d)",
                    Breakend.gene(), SvData.id(), SvData.type(),
                    SvData.startChromosome(), SvData.startPosition(), SvData.endChromosome(), SvData.endPosition(),
                    mComparisonStartGenomePosition.chromosome(), mComparisonStartGenomePosition.position(),
                    mComparisonEndGenomePosition.chromosome(), mComparisonEndGenomePosition.position());
        }
        else
        {
            return String.format("%s %d_%s %s:%d-%s:%d",
                    Breakend.gene(), SvData.id(), SvData.type(),
                    SvData.startChromosome(), SvData.startPosition(), SvData.endChromosome(), SvData.endPosition());
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("%s", Breakend.reportedDisruption()));
        values.add(String.format("%s", Breakend.regionType()));
        values.add(String.format("%s", Breakend.codingType()));
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

        if(!otherSv.SvData.startChromosome().equals(mComparisonStartGenomePosition.chromosome()) || !otherSv.SvData.endChromosome().equals(mComparisonEndGenomePosition.chromosome()))
            return false;

        if(otherSv.SvData.startPosition() != mComparisonStartGenomePosition.position() || otherSv.SvData.endPosition() != mComparisonEndGenomePosition.position())
            return false;

        if(otherSv.SvData.startOrientation() != SvData.startOrientation() || otherSv.SvData.endOrientation() != SvData.endOrientation())
            return false;

        if(!otherSv.Breakend.gene().equals(Breakend.gene()))
            return false;

        if((otherSv.mCheckTranscript || mCheckTranscript) && !otherSv.Breakend.transcriptId().equals(Breakend.transcriptId()))
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final DisruptionData otherBreakend = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_REGION_TYPE, Breakend.regionType().toString(), otherBreakend.Breakend.regionType().toString());
        checkDiff(diffs, FLD_CODING_CONTEXT, Breakend.codingType().toString(), otherBreakend.Breakend.codingType().toString());
        checkDiff(diffs, FLD_REPORTED, reportable(), otherBreakend.reportable());
        checkDiff(diffs, FLD_GENE_ORIENT, Breakend.geneOrientation(), otherBreakend.Breakend.geneOrientation());
        checkDiff(diffs, FLD_NEXT_SPLICE, Breakend.nextSpliceExonRank(), otherBreakend.Breakend.nextSpliceExonRank());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
