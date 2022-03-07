package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.CommonUtils.filtersStr;
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
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;
import com.hartwig.hmftools.compar.MismatchType;

public class SomaticVariantData implements ComparableItem
{
    public final SomaticVariant Variant;
    public final Set<String> Filters;

    public SomaticVariantData(final SomaticVariant variant)
    {
        Variant = variant;
        Filters = Arrays.stream(variant.filter().split(";", -1)).collect(Collectors.toSet());
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public String key()
    {
        return String.format("%s:%d %s>%s %s",
                Variant.chromosome(), Variant.position(), Variant.ref(), Variant.alt(), Variant.type());
    }

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
        final SomaticVariantData otherVar = (SomaticVariantData)other;

        if(!Variant.chromosome().equals(otherVar.Variant.chromosome()) || Variant.position() != otherVar.Variant.position())
            return false;

        if(!Variant.ref().equals(otherVar.Variant.ref()) || !Variant.alt().equals(otherVar.Variant.alt()))
            return false;

        if(Variant.type() != otherVar.Variant.type())
            return false;

        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        final SomaticVariantData otherVar = (SomaticVariantData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Variant.reported(), otherVar.Variant.reported());

        if(matchLevel != REPORTABLE)
        {
            checkDiff(diffs, "hotspot", Variant.hotspot().toString(), otherVar.Variant.hotspot().toString());
            checkDiff(diffs, "tier", Variant.tier().toString(), otherVar.Variant.tier().toString());
            checkDiff(diffs, "biallelic", Variant.biallelic(), otherVar.Variant.biallelic());
            checkDiff(diffs, "gene", Variant.gene(), otherVar.Variant.gene());
            checkDiff(diffs, "canonicalEffect", Variant.canonicalEffect(), otherVar.Variant.canonicalEffect());
            checkDiff(diffs, "canonicalCodingEffect", Variant.canonicalCodingEffect().toString(), otherVar.Variant.canonicalCodingEffect().toString());
            checkDiff(diffs, "canonicalHgvsCoding", Variant.canonicalHgvsCodingImpact(), otherVar.Variant.canonicalHgvsCodingImpact());
            checkDiff(diffs, "canonicalHgvsProtein", Variant.canonicalHgvsProteinImpact(), otherVar.Variant.canonicalHgvsProteinImpact());
            checkDiff(diffs, "otherReportedEffects", Variant.otherReportedEffects(), otherVar.Variant.otherReportedEffects());
            checkDiff(diffs, "subclonalLikelihood", Variant.subclonalLikelihood(), otherVar.Variant.subclonalLikelihood());
            checkDiff(diffs, "qual", (int)Variant.qual(), (int)otherVar.Variant.qual());
        }

        if(diffs.isEmpty())
            return null;

        return new Mismatch(this, other, VALUE, diffs);
    }

    @Override
    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final SomaticVariantData otherVar = (SomaticVariantData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Variant.reported(), otherVar.Variant.reported());

        if(matchLevel == REPORTABLE)
            return diffs;

        // key fields
        checkDiff(diffs, "hotspot", Variant.hotspot().toString(), otherVar.Variant.hotspot().toString());
        checkDiff(diffs, "tier", Variant.tier().toString(), otherVar.Variant.tier().toString());
        checkDiff(diffs, "biallelic", Variant.biallelic(), otherVar.Variant.biallelic());
        checkDiff(diffs, "gene", Variant.gene(), otherVar.Variant.gene());
        checkDiff(diffs, "canonicalEffect", Variant.canonicalEffect(), otherVar.Variant.canonicalEffect());
        checkDiff(diffs, "canonicalCodingEffect", Variant.canonicalCodingEffect().toString(), otherVar.Variant.canonicalCodingEffect().toString());
        checkDiff(diffs, "canonicalHgvsCoding", Variant.canonicalHgvsCodingImpact(), otherVar.Variant.canonicalHgvsCodingImpact());
        checkDiff(diffs, "canonicalHgvsProtein", Variant.canonicalHgvsProteinImpact(), otherVar.Variant.canonicalHgvsProteinImpact());
        checkDiff(diffs, "otherReportedEffects", Variant.otherReportedEffects(), otherVar.Variant.otherReportedEffects());
        checkDiff(diffs, "subclonalLikelihood", Variant.subclonalLikelihood(), otherVar.Variant.subclonalLikelihood());
        checkDiff(diffs, "qual", (int)Variant.qual(), (int)otherVar.Variant.qual());

        // compare filters
        Set<String> origFilters = Filters;
        Set<String> newFilters = otherVar.Filters;

        Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
        Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

        if(!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty())
        {
            diffs.add(String.format("%s(%.3f/%.3f)", "filter", filtersStr(origFilterDiffs), filtersStr(newFilterDiffs)));
        }

        return diffs;
    }

}
