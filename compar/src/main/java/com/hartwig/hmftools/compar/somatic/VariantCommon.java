package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_LPS;
import static com.hartwig.hmftools.compar.somatic.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.compar.DiffThresholds;

public final class VariantCommon
{
    protected static final String FLD_HOTSPOT = "Hotspot";
    protected static final String FLD_TIER = "Tier";
    protected static final String FLD_BIALLELIC = "Biallelic";
    protected static final String FLD_GENE = "Gene";
    protected static final String FLD_CANON_EFFECT = "CanonicalEffect";
    protected static final String FLD_CODING_EFFECT = "CanonicalCodingEffect";
    protected static final String FLD_HGVS_CODING = "CanonicalHgvsCoding";
    protected static final String FLD_HGVS_PROTEIN = "CanonicalHgvsProtein";
    protected static final String FLD_OTHER_REPORTED = "OtherReportedEffects";
    protected static final String FLD_QUAL = "Qual";

    protected static final List<String> findVariantDiffs(
            final SomaticVariant refVar, final SomaticVariant otherVar, final DiffThresholds thresholds)
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

    protected static List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_BIALLELIC, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL);
    }

    protected static void addDisplayValues(final SomaticVariant variant, final List<String> values)
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
