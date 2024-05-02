package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_SV;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineSvData implements ComparableItem
{
    public final LinxGermlineSv SvData;
    private final boolean mIsReported;
    private final GenomePositionImpl mComparisonStartGenomePosition;
    private final GenomePositionImpl mComparisonEndGenomePosition;

    protected static final String FLD_GERMLINE_FRAGS = "GermlineFragments";

    public GermlineSvData(final LinxGermlineSv svData, boolean isReported, final GenomePositionImpl comparisonStartGenomePosition,
            final GenomePositionImpl comparisonEndGenomePosition)
    {
        SvData = svData;
        mIsReported = isReported;
        mComparisonStartGenomePosition = comparisonStartGenomePosition;
        mComparisonEndGenomePosition = comparisonEndGenomePosition;
    }

    @Override
    public Category category() { return GERMLINE_SV; }

    @Override
    public String key()
    {
        if(mComparisonStartGenomePosition.position() != SvData.PositionStart || mComparisonEndGenomePosition.position() != SvData.PositionEnd)
        {
            return String.format("%s:%s %s:%d:%d-%s:%d%d %s liftover(%s:%d-%s:%d)",
                    SvData.EventId, SvData.Type, SvData.ChromosomeStart, SvData.PositionStart, SvData.OrientStart,
                    SvData.ChromosomeEnd, SvData.PositionEnd, SvData.OrientEnd, SvData.GeneName,
                    mComparisonStartGenomePosition.chromosome(), mComparisonStartGenomePosition.position(),
                    mComparisonEndGenomePosition.chromosome(), mComparisonEndGenomePosition.position());
        }
        else
        {
            return String.format("%s:%s %s:%d:%d-%s:%d%d %s",
                    SvData.EventId, SvData.Type, SvData.ChromosomeStart, SvData.PositionStart, SvData.OrientStart,
                    SvData.ChromosomeEnd, SvData.PositionEnd, SvData.OrientEnd, SvData.GeneName);
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("%s", mIsReported));
        values.add(String.format("%d", SvData.GermlineFragments));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return mIsReported;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineSvData otherSv = (GermlineSvData)other;

        if(otherSv.SvData.Type != SvData.Type)
            return false;

        if(!otherSv.SvData.ChromosomeStart.equals(mComparisonStartGenomePosition.chromosome()) || !otherSv.SvData.ChromosomeEnd.equals(mComparisonEndGenomePosition.chromosome()))
            return false;

        if(otherSv.SvData.PositionStart != mComparisonStartGenomePosition.position() || otherSv.SvData.PositionEnd != mComparisonEndGenomePosition.position())
            return false;

        if(otherSv.SvData.OrientStart != SvData.OrientStart || otherSv.SvData.OrientEnd != SvData.OrientEnd)
            return false;

        if(!otherSv.SvData.GeneName.equals(SvData.GeneName))
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final GermlineSvData otherSv = (GermlineSvData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_QUAL, (int) SvData.QualScore, (int) otherSv.SvData.QualScore, thresholds);
        checkDiff(diffs, FLD_GERMLINE_FRAGS, SvData.GermlineFragments, otherSv.SvData.GermlineFragments, thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
