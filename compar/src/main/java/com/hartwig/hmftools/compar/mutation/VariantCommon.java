package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;

import java.util.List;

import com.google.common.collect.Lists;

public final class VariantCommon
{
    static final String FLD_HOTSPOT = "Hotspot";
    static final String FLD_TIER = "Tier";
    static final String FLD_BIALLELIC = "Biallelic";
    static final String FLD_GENE = "Gene";
    static final String FLD_CANON_EFFECT = "CanonicalEffect";
    static final String FLD_CODING_EFFECT = "CanonicalCodingEffect";
    static final String FLD_HGVS_CODING = "CanonicalHgvsCoding";
    static final String FLD_HGVS_PROTEIN = "CanonicalHgvsProtein";
    static final String FLD_OTHER_REPORTED = "OtherReportedEffects";
    static final String FLD_VARIANT_COPY_NUMBER = "VariantCopyNumber";
    static final String FLD_PURITY_ADJUSTED_VAF = "PurityAdjustedVaf";
    static final String FLD_TUMOR_SUPPORTING_READ_COUNT = "TumorSupportingReadCount";
    static final String FLD_TUMOR_TOTAL_READ_COUNT = "TumorTotalReadCount";

    static List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_BIALLELIC, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL, FLD_VARIANT_COPY_NUMBER, FLD_PURITY_ADJUSTED_VAF,
                FLD_TUMOR_SUPPORTING_READ_COUNT, FLD_TUMOR_TOTAL_READ_COUNT);
    }
}
