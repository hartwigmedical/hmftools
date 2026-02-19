package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record CopyNumberData(
        String chromosome,
        int positionStart,
        int positionEnd,
        double copyNumber,
        double majorAlleleCopyNumber,
        CopyNumberMethod method,
        BasePosition comparisonPositionStart,
        BasePosition comparisonPositionEnd
) implements ComparableItem
{
    static final String FLD_COPY_NUMBER = "CopyNumber";
    static final String FLD_MAJOR_ALLELE_CN = "MajorAlleleCopyNumber";
    static final String FLD_METHOD = "Method";

    public CategoryType category() {
        return COPY_NUMBER;
    }

    @Override
    public String key()
    {
        if(comparisonPositionStart.equals(new BasePosition(chromosome, positionStart))
                && comparisonPositionEnd.equals(new BasePosition(chromosome, positionEnd)))
        {
            return format("%s:%d_%d", chromosome, positionStart, positionEnd);
        }
        else
        {
            return format("%s:%d_%d liftover(%s_%s)", chromosome, positionStart, positionEnd,
                    comparisonPositionStart, comparisonPositionEnd);
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", copyNumber));
        values.add(format("%.2f", majorAlleleCopyNumber));
        values.add(format("%s", method));
        return values;
    }

    @Override
    public boolean reportable() {
        return false;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        if(!comparisonPositionStart.Chromosome.equals(otherCn.chromosome) || !comparisonPositionEnd.Chromosome.equals(otherCn.chromosome))
            return false;

        return comparisonPositionStart.Position == otherCn.positionStart && comparisonPositionEnd.Position == otherCn.positionEnd;
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_COPY_NUMBER, copyNumber, otherCn.copyNumber, thresholds);
        checkDiff(diffs, FLD_MAJOR_ALLELE_CN, majorAlleleCopyNumber, otherCn.majorAlleleCopyNumber, thresholds);
        checkDiff(diffs, FLD_METHOD, method.toString(), otherCn.method.toString());

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
