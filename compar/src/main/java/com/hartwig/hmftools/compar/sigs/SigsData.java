package com.hartwig.hmftools.compar.sigs;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sigs.SignatureAllocationFile.PERCENT_FLD;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record SigsData(SignatureAllocation SignatureAllocation) implements ComparableItem
{
    public static String FLD_PERCENT = capitalize(PERCENT_FLD);

    @Override
    public CategoryType category()
    {
        return CategoryType.SIGS;
    }

    @Override
    public String key()
    {
        return SignatureAllocation.signature();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", SignatureAllocation.percent()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return false;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SigsData otherData = (SigsData)other;

        return otherData.SignatureAllocation.signature().equals(SignatureAllocation.signature());
    }


    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final SignatureAllocation otherData = ((SigsData) other).SignatureAllocation;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_PERCENT, SignatureAllocation.percent(), otherData.percent(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
