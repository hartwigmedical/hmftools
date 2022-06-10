package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.DiffFunctions.checkFilterDiffs;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.findVariantDiffs;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.variantsMatch;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class GermlineVariantData implements ComparableItem
{
    public final SomaticVariant Variant;
    public final Set<String> Filters;

    public GermlineVariantData(final SomaticVariant variant)
    {
        Variant = variant;
        Filters = Arrays.stream(variant.filter().split(";", -1)).collect(Collectors.toSet());
    }

    @Override
    public Category category() { return GERMLINE_VARIANT; }

    @Override
    public String key() { return SomaticVariantData.key(Variant); }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("Qual(%.0f)", Variant.qual()));
        values.add(String.format("Tier(%s)", Variant.tier().toString()));
        values.add(String.format("TotalReadCount(%d)", Variant.totalReadCount()));
        values.add(String.format("AlleleReadCount(%d)", Variant.alleleReadCount()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return Variant.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineVariantData otherVar = (GermlineVariantData) other;
        return variantsMatch(Variant, otherVar.Variant);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final GermlineVariantData otherVar = (GermlineVariantData) other;

        final List<String> diffs = findVariantDiffs(Variant, otherVar.Variant, thresholds);

        checkFilterDiffs(Filters, otherVar.Filters, diffs);

        if(diffs.isEmpty())
            return null;

        return new Mismatch(this, other, VALUE, diffs);
    }

}
