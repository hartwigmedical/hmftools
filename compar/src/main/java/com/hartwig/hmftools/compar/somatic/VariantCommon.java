package com.hartwig.hmftools.compar.somatic;

import static com.hartwig.hmftools.compar.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;

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

    protected static List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_BIALLELIC, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL);
    }
}
