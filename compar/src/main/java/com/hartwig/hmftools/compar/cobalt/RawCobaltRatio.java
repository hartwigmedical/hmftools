package com.hartwig.hmftools.compar.cobalt;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.cobalt.CobaltRatioComparer.REFERENCE_GC_DIPLOID_RATIO;
import static com.hartwig.hmftools.compar.cobalt.CobaltRatioComparer.REFERENCE_GC_RATIO;
import static com.hartwig.hmftools.compar.cobalt.CobaltRatioComparer.REFERENCE_READ_COUNT;
import static com.hartwig.hmftools.compar.cobalt.CobaltRatioComparer.TUMOR_GC_RATIO;
import static com.hartwig.hmftools.compar.cobalt.CobaltRatioComparer.TUMOR_READ_COUNT;
import static com.hartwig.hmftools.compar.common.Category.COBALT_RATIO;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

import org.jetbrains.annotations.NotNull;

public record RawCobaltRatio(
        @NotNull String chromosome,
        int position,
        double referenceReadCount,
        double tumorReadCount,
        double referenceGcRatio,
        double tumorGcRatio,
        double referenceGcDiploidRatio
) implements ComparableItem
{
    @Override
    public Category category()
    {
        return COBALT_RATIO;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        if(!(other instanceof final RawCobaltRatio otherRC))
        {
            return false;
        }
        return chromosome.equals(otherRC.chromosome) && position == otherRC.position;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final RawCobaltRatio otherCn = (RawCobaltRatio) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, REFERENCE_READ_COUNT, referenceReadCount, otherCn.referenceReadCount, thresholds);
        checkDiff(diffs, TUMOR_READ_COUNT, tumorReadCount, otherCn.tumorReadCount, thresholds);
        checkDiff(diffs, REFERENCE_GC_RATIO, referenceGcRatio, otherCn.referenceGcRatio, thresholds);
        checkDiff(diffs, TUMOR_GC_RATIO, tumorGcRatio, otherCn.tumorGcRatio, thresholds);
        checkDiff(diffs, REFERENCE_GC_DIPLOID_RATIO, referenceGcDiploidRatio, otherCn.referenceGcDiploidRatio, thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    @Override
    public String key()
    {
        return format("%s:%d", chromosome, position);
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.4f", referenceReadCount));
        values.add(format("%.4f", tumorReadCount));
        values.add(format("%.4f", referenceGcRatio));
        values.add(format("%.4f", tumorGcRatio));
        values.add(format("%.4f", referenceGcDiploidRatio));
        return values;
    }
}
