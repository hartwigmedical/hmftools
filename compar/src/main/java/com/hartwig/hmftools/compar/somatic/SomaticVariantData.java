package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.DiffFunctions.checkFilterDiffs;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

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

public class SomaticVariantData implements ComparableItem
{
    public final SomaticVariant Variant;
    public final Set<String> Filters;

    protected static final String FLD_QUAL = "qual";
    protected static final String FLD_SUBCLONAL_LIKELIHOOD = "subclonalLikelihood";

    public SomaticVariantData(final SomaticVariant variant)
    {
        Variant = variant;
        Filters = Arrays.stream(variant.filter().split(";", -1)).collect(Collectors.toSet());
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public String key() { return key(Variant); }

    protected static String key(final SomaticVariant var)
    {
        return String.format("%s:%d %s>%s %s", var.chromosome(), var.position(), var.ref(), var.alt(), var.type());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("Qual(%.0f)", Variant.qual()));
        values.add(String.format("Tier(%s)", Variant.tier().toString()));
        values.add(String.format("TotalReadCount(%d)", Variant.totalReadCount()));
        values.add(String.format("AlleleReadCount(%d)", Variant.alleleReadCount()));
        values.add(String.format("LPS(%s)", Variant.localPhaseSetsStr()));
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
        final SomaticVariantData otherVar = (SomaticVariantData) other;
        return variantsMatch(Variant, otherVar.Variant);
    }

    protected static boolean variantsMatch(final SomaticVariant var1, final SomaticVariant var2)
    {
        if(!var1.chromosome().equals(var1.chromosome()) || var1.position() != var2.position())
            return false;

        if(!var1.ref().equals(var2.ref()) || !var1.alt().equals(var2.alt()))
            return false;

        if(var1.type() != var2.type())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final SomaticVariantData otherVar = (SomaticVariantData) other;

        final List<String> diffs = findVariantDiffs(Variant, otherVar.Variant, thresholds);

        checkDiff(diffs, "hasLPS", Variant.hasLocalPhaseSets(), otherVar.Variant.hasLocalPhaseSets());

        checkDiff(diffs, FLD_SUBCLONAL_LIKELIHOOD, Variant.subclonalLikelihood(), otherVar.Variant.subclonalLikelihood(), thresholds);

        // compare filters
        checkFilterDiffs(Filters, otherVar.Filters, diffs);

        if(diffs.isEmpty())
            return null;

        return new Mismatch(this, other, VALUE, diffs);
    }

    // shared with germline variant
    public static final List<String> findVariantDiffs(
            final SomaticVariant refVar, final SomaticVariant otherVar, final DiffThresholds thresholds)
    {
        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", refVar.reported(), otherVar.reported());

        checkDiff(diffs, "hotspot", refVar.hotspot().toString(), otherVar.hotspot().toString());
        checkDiff(diffs, "tier", refVar.tier().toString(), otherVar.tier().toString());
        checkDiff(diffs, "biallelic", refVar.biallelic(), otherVar.biallelic());
        checkDiff(diffs, "gene", refVar.gene(), otherVar.gene());
        checkDiff(diffs, "canonicalEffect", refVar.canonicalEffect(), otherVar.canonicalEffect());
        checkDiff(diffs, "canonicalCodingEffect", refVar.canonicalCodingEffect().toString(), otherVar.canonicalCodingEffect()
                .toString());
        checkDiff(diffs, "canonicalHgvsCoding", refVar.canonicalHgvsCodingImpact(), otherVar.canonicalHgvsCodingImpact());
        checkDiff(diffs, "canonicalHgvsProtein", refVar.canonicalHgvsProteinImpact(), otherVar.canonicalHgvsProteinImpact());
        checkDiff(diffs, "otherReportedEffects", refVar.otherReportedEffects(), otherVar.otherReportedEffects());

        checkDiff(diffs, FLD_QUAL, (int) refVar.qual(), (int) otherVar.qual(), thresholds);

        return diffs;
    }
}
