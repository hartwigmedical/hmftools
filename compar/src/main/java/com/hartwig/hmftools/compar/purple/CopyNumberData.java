package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.COPY_NUMBER;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

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
}
