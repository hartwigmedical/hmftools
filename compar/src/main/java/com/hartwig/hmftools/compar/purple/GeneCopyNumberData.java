package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GeneCopyNumberData implements ComparableItem
{
    public final GeneCopyNumber CopyNumber;

    protected static final String FLD_MIN_COPY_NUMBER = "MinCopyNumber";
    protected static final String FLD_MAX_COPY_NUMBER = "MaxCopyNumber";
    protected static final String FLD_MIN_REGION_START = "MinRegionStart";
    protected static final String FLD_MIN_REGION_END = "MinRegionEnd";

    public GeneCopyNumberData(final GeneCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    public CategoryType category() {
        return GENE_COPY_NUMBER;
    }

    @Override
    public String key()
    {
        return format("%s", CopyNumber.geneName());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", CopyNumber.minCopyNumber()));
        values.add(format("%.2f", CopyNumber.maxCopyNumber()));
        values.add(format("%d", CopyNumber.MinRegionStart));
        values.add(format("%d", CopyNumber.MinRegionEnd));
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
    public String geneName() { return CopyNumber.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GeneCopyNumberData otherCn = (GeneCopyNumberData)other;

        return CopyNumber.geneName().equals(otherCn.CopyNumber.geneName());
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches)
    {
        final GeneCopyNumberData otherCn = (GeneCopyNumberData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_MIN_COPY_NUMBER, CopyNumber.minCopyNumber(), otherCn.CopyNumber.minCopyNumber(), thresholds);
        checkDiff(diffs, FLD_MAX_COPY_NUMBER, CopyNumber.maxCopyNumber(), otherCn.CopyNumber.maxCopyNumber(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
