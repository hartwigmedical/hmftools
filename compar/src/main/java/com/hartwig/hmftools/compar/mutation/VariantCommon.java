package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_BIALLELIC_PROB;
import static com.hartwig.hmftools.compar.mutation.SomaticVariantData.FLD_SUBCLONAL_LIKELIHOOD;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;

public final class VariantCommon
{
    static final String FLD_HOTSPOT = "Hotspot";
    static final String FLD_TIER = "Tier";
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

    public static void registerThresholds(final CategoryType category, final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(category, FLD_QUAL, 50, 0.2);
        thresholds.addFieldThreshold(category, FLD_SUBCLONAL_LIKELIHOOD, 0.6, 0);
        thresholds.addFieldThreshold(category, FLD_BIALLELIC_PROB, 0.3, 0);
        thresholds.addFieldThreshold(category, FLD_VARIANT_COPY_NUMBER, 0.3, 0.3);
        thresholds.addFieldThreshold(category, FLD_PURITY_ADJUSTED_VAF, 0.2, 0);
        thresholds.addFieldThreshold(category, FLD_TUMOR_SUPPORTING_READ_COUNT, 1, 0.2);
        thresholds.addFieldThreshold(category, FLD_TUMOR_TOTAL_READ_COUNT, 1, 0.2);
    }

    static List<String> comparedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL, FLD_VARIANT_COPY_NUMBER, FLD_PURITY_ADJUSTED_VAF,
                FLD_TUMOR_SUPPORTING_READ_COUNT, FLD_TUMOR_TOTAL_READ_COUNT);
    }
}
