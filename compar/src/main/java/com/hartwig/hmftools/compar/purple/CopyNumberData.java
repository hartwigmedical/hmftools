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

public class CopyNumberData implements ComparableItem
{
    public final PurpleCopyNumber CopyNumber;

    public CopyNumberData(final PurpleCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    public Category category() { return COPY_NUMBER; }

    public boolean reportable() { return false; }

    public boolean matches(final ComparableItem other)
    {
        final CopyNumberData otherCn = (CopyNumberData)other;

        if(CopyNumber.chromosome() != otherCn.CopyNumber.chromosome())
            return false;

        if(CopyNumber.start() != otherCn.CopyNumber.start() || CopyNumber.end() != otherCn.CopyNumber.end())
            return false;

        return true;
    }

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

    public String description()
    {
        return String.format("location(%s %d-%d) %s_%s",
                CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end(),
                CopyNumber.segmentStartSupport(), CopyNumber.segmentEndSupport());
    }
}
