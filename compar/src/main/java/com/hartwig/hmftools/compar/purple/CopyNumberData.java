package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class CopyNumberData implements ComparableItem
{
    public final PurpleCopyNumber CopyNumber;

    protected static final String FLD_COPY_NUMBER = "CopyNumber";
    protected static final String FLD_MAJOR_ALLELE_CN = "MajorAlleleCopyNumber";
    protected static final String FLD_METHOD = "Method";

    public CopyNumberData(final PurpleCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    public Category category() {
        return COPY_NUMBER;
    }

    @Override
    public String key()
    {
        return format("%s:%d_%d", CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end());
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
    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData)other;

        return CopyNumber.chromosome().equals(otherCn.CopyNumber.chromosome())
                && CopyNumber.start() == otherCn.CopyNumber.start() && CopyNumber.end() == otherCn.CopyNumber.end();
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final CopyNumberData otherCn = (CopyNumberData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_COPY_NUMBER, CopyNumber.averageTumorCopyNumber(), otherCn.CopyNumber.averageTumorCopyNumber(), thresholds);
        checkDiff(diffs, FLD_COPY_NUMBER, CopyNumber.majorAlleleCopyNumber(), otherCn.CopyNumber.majorAlleleCopyNumber(), thresholds);
        checkDiff(diffs, FLD_COPY_NUMBER, CopyNumber.method().toString(), otherCn.CopyNumber.method().toString());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
