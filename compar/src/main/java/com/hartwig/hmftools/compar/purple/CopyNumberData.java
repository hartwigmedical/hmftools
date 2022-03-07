package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.COPY_NUMBER;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class CopyNumberData implements ComparableItem
{
    public final PurpleCopyNumber CopyNumber;

    public CopyNumberData(final PurpleCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    @Override
    public Category category() { return COPY_NUMBER; }

    @Override
    public String key()
    {
        return String.format("location(%s %d-%d) %s_%s",
                CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end(),
                CopyNumber.segmentStartSupport(), CopyNumber.segmentEndSupport());
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
    public boolean reportable() { return false; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData)other;

        if(CopyNumber.chromosome() != otherCn.CopyNumber.chromosome())
            return false;

        if(CopyNumber.start() != otherCn.CopyNumber.start() || CopyNumber.end() != otherCn.CopyNumber.end())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        return null;
    }

    @Override
    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final CopyNumberData otherCn = (CopyNumberData)other;

        final List<String> diffs = Lists.newArrayList();

        if(matchLevel == REPORTABLE)
            return diffs;

        checkDiff(diffs, "segStartSupport", CopyNumber.segmentStartSupport().toString(), otherCn.CopyNumber.segmentStartSupport().toString());
        checkDiff(diffs, "segEndSupport", CopyNumber.segmentEndSupport().toString(), otherCn.CopyNumber.segmentEndSupport().toString());
        checkDiff(diffs, "method", CopyNumber.method().toString(), otherCn.CopyNumber.method().toString());

        return diffs;
    }
}
