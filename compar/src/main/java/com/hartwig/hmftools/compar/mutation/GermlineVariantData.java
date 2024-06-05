package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.Category.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkFilterDiffs;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_BIALLELIC;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CANON_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CODING_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_GENE;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_CODING;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_PROTEIN;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HOTSPOT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_OTHER_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TIER;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineVariantData implements ComparableItem
{
    public final GermlineVariant Variant;
    public final Set<String> Filters;
    private final BasePosition mComparisonPosition;

    public GermlineVariantData(final GermlineVariant variant, final BasePosition comparisonPosition)
    {
        Variant = variant;
        Filters = Arrays.stream(variant.filter().split(";", -1)).collect(Collectors.toSet());
        mComparisonPosition = comparisonPosition;
    }

    @Override
    public Category category() { return GERMLINE_VARIANT; }

    @Override
    public String key()
    {
        if(mComparisonPosition.Position != Variant.position())
        {
            return String.format("%s:%d %s>%s %s liftover(%s)",
                    Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type(), mComparisonPosition);
        }
        else
        {
            return String.format("%s:%d %s>%s %s", Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type());
        }
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        addDisplayValues(Variant, values);
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

        if(!mComparisonPosition.Chromosome.equals(otherVar.Variant.chromosome()) || mComparisonPosition.Position != otherVar.Variant.position())
            return false;

        if(!Variant.ref().equals(otherVar.Variant.ref()) || !Variant.alt().equals(otherVar.Variant.alt()))
            return false;

        if(Variant.type() != otherVar.Variant.type())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final GermlineVariantData otherVar = (GermlineVariantData) other;

        final List<String> diffs = findVariantDiffs(Variant, otherVar.Variant, thresholds);

        checkFilterDiffs(Filters, otherVar.Filters, diffs);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    private static List<String> findVariantDiffs(
            final GermlineVariant refVar, final GermlineVariant otherVar, final DiffThresholds thresholds)
    {
        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_REPORTED, refVar.reported(), otherVar.reported());
        checkDiff(diffs, FLD_HOTSPOT, refVar.hotspot().toString(), otherVar.hotspot().toString());
        checkDiff(diffs, FLD_TIER, refVar.tier().toString(), otherVar.tier().toString());
        checkDiff(diffs, FLD_BIALLELIC, refVar.biallelic(), otherVar.biallelic());
        checkDiff(diffs, FLD_GENE, refVar.gene(), otherVar.gene());
        checkDiff(diffs, FLD_CANON_EFFECT, refVar.canonicalEffect(), otherVar.canonicalEffect());
        checkDiff(diffs, FLD_CODING_EFFECT, refVar.canonicalCodingEffect().toString(), otherVar.canonicalCodingEffect()
                .toString());
        checkDiff(diffs, FLD_HGVS_CODING, refVar.canonicalHgvsCodingImpact(), otherVar.canonicalHgvsCodingImpact());
        checkDiff(diffs, FLD_HGVS_PROTEIN, refVar.canonicalHgvsProteinImpact(), otherVar.canonicalHgvsProteinImpact());
        checkDiff(diffs, FLD_OTHER_REPORTED, refVar.otherReportedEffects(), otherVar.otherReportedEffects());

        checkDiff(diffs, FLD_QUAL, (int) refVar.qual(), (int) otherVar.qual(), thresholds);

        return diffs;
    }

    protected static void addDisplayValues(final GermlineVariant variant, final List<String> values)
    {
        values.add(String.format("%s", variant.reported()));
        values.add(String.format("%s", variant.hotspot()));
        values.add(String.format("%s", variant.tier()));
        values.add(String.format("%s", variant.biallelic()));
        values.add(String.format("%s", variant.gene()));
        values.add(String.format("%s", variant.canonicalEffect()));
        values.add(String.format("%s", variant.canonicalCodingEffect()));
        values.add(String.format("%s", variant.canonicalHgvsCodingImpact()));
        values.add(String.format("%s", variant.canonicalHgvsProteinImpact()));
        values.add(String.format("%s", variant.otherReportedEffects()));
        values.add(String.format("%.0f", variant.qual()));
    }
}
