package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.compar.DiffThresholds;

public final class VariantCommon
{
    protected static final String FLD_QUAL = "qual";

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

    public static void addDisplayValues(final SomaticVariant variant, final List<String> values)
    {
        values.add(String.format("qual=%.0f", variant.qual()));
        values.add(String.format("hotspot=%s", variant.hotspot()));
        values.add(String.format("tier=%s", variant.tier()));
        values.add(String.format("biallelic=%s", variant.biallelic()));
        values.add(String.format("gene=%s", variant.gene()));
        values.add(String.format("canonicalEffect=%s", variant.canonicalEffect()));
        values.add(String.format("canonicalCodingEffect=%s", variant.canonicalCodingEffect()));
        values.add(String.format("canonicalHgvsCoding=%s", variant.canonicalHgvsCodingImpact()));
        values.add(String.format("canonicalHgvsProtein=%s", variant.canonicalHgvsProteinImpact()));
        values.add(String.format("hasLPS=%s", variant.hasLocalPhaseSets()));
    }

}
