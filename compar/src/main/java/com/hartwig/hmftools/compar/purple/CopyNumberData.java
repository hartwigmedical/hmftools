package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.COPY_NUMBER;

import java.util.List;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class CopyNumberData extends ComparableItem
{
    public final String Chromosome;
    public int PositionStart;
    public int PositionEnd;
    public double CopyNumber;
    public double MajorAlleleCopyNumber;
    final CopyNumberMethod Method;

    public CopyNumberData(
            final String chromosome, int positionStart, int positionEnd, double copyNumber, double majorAlleleCopyNumber,
            final CopyNumberMethod method, final List<FieldInfo> fields)
    {
        Chromosome = chromosome;
        PositionStart = positionStart;
        PositionEnd = positionEnd;
        CopyNumber = copyNumber;
        Method = method;

        addDoubleValue(CopyNumberComparer.Fields.CopyNumber.toString(), copyNumber, fields);
        addDoubleValue(CopyNumberComparer.Fields.MajorAlleleCopyNumber.toString(), majorAlleleCopyNumber, fields);
        addStringValue(CopyNumberComparer.Fields.Method.toString(), method.toString(), fields);
    }

    public CategoryType category() {
        return COPY_NUMBER;
    }

    @Override
    public String key()
    {
        return format("%s:%d_%d", Chromosome, PositionStart, PositionEnd);
    }

    @Override
    public boolean reportable() {
        return false;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        return Chromosome.equals(otherCn.Chromosome) && PositionStart == otherCn.PositionStart && PositionEnd == otherCn.PositionEnd;
    }
}
