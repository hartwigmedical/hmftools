package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class CopyNumberData implements ComparableItem
{
    public final PurpleCopyNumber CopyNumber;
    public final BasePosition mComparisonPositionStart;
    public final BasePosition mComparisonPositionEnd;

    protected static final String FLD_COPY_NUMBER = "CopyNumber";
    protected static final String FLD_MAJOR_ALLELE_CN = "MajorAlleleCopyNumber";
    protected static final String FLD_METHOD = "Method";

    public CopyNumberData(final PurpleCopyNumber copyNumber, final BasePosition comparisonPositionStart,
            final BasePosition comparisonPositionEnd)
    {
        CopyNumber = copyNumber;
        mComparisonPositionStart = comparisonPositionStart;
        mComparisonPositionEnd = comparisonPositionEnd;
    }

    public Category category() {
        return COPY_NUMBER;
    }

    @Override
    public String key()
    {
        if(mComparisonPositionStart.equals(new BasePosition(CopyNumber.chromosome(), CopyNumber.start()))
                && mComparisonPositionEnd.equals(new BasePosition(CopyNumber.chromosome(), CopyNumber.end())))
        {
            return format("%s:%d_%d", CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end());
        }
        else
        {
            return format("%s:%d_%d liftover(%s_%s)", CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end(),
                    mComparisonPositionStart, mComparisonPositionEnd);
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", CopyNumber.averageTumorCopyNumber()));
        values.add(format("%.2f", CopyNumber.majorAlleleCopyNumber()));
        values.add(format("%s", CopyNumber.method()));
        return values;
    }

    @Override
    public boolean reportable() {
        return false;
    }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        if(!mComparisonPositionStart.Chromosome.equals(otherCn.CopyNumber.chromosome())
                || !mComparisonPositionEnd.Chromosome.equals(otherCn.CopyNumber.chromosome()))
            return false;

        return mComparisonPositionStart.Position == otherCn.CopyNumber.start()
                && mComparisonPositionEnd.Position == otherCn.CopyNumber.end();
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_COPY_NUMBER, CopyNumber.averageTumorCopyNumber(), otherCn.CopyNumber.averageTumorCopyNumber(), thresholds);
        checkDiff(diffs, FLD_MAJOR_ALLELE_CN, CopyNumber.majorAlleleCopyNumber(), otherCn.CopyNumber.majorAlleleCopyNumber(), thresholds);
        checkDiff(diffs, FLD_METHOD, CopyNumber.method().toString(), otherCn.CopyNumber.method().toString());

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
