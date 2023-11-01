package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;

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
