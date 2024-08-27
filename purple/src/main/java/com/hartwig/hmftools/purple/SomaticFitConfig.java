package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleConstants.HIGHLY_DIPLOID_PERCENTAGE_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_MIN_PEAK_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_MIN_PURITY_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_MIN_PURITY_SPREAD_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_MIN_VARIANTS_DEFAULT;
import static com.hartwig.hmftools.purple.PurpleConstants.SOMATIC_PENALTY_WEIGHT_DEFAULT;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SomaticFitConfig
{
    public static final String SOMATIC_MIN_PEAK = "somatic_min_peak";
    public static final String SOMATIC_MIN_TOTAL = "somatic_min_variants";
    public static final String SOMATIC_MIN_PURITY = "somatic_min_purity";
    public static final String SOMATIC_MIN_PURITY_SPREAD = "somatic_min_purity_spread";
    public static final String SOMATIC_PENALTY_WEIGHT = "somatic_penalty_weight";
    public static final String HIGHLY_DIPLOID_PERCENTAGE = "highly_diploid_percentage";

    public int MinTotalVariants;
    public int MinPeakVariants;
    public double MinSomaticPurity;
    public double MinSomaticPuritySpread;
    public double SomaticPenaltyWeight;
    public double HighlyDiploidPercentage;

    public SomaticFitConfig(final ConfigBuilder configBuilder)
    {
        MinTotalVariants = configBuilder.getInteger(SOMATIC_MIN_TOTAL);
        MinPeakVariants = configBuilder.getInteger(SOMATIC_MIN_PEAK);
        MinSomaticPurity = configBuilder.getDecimal(SOMATIC_MIN_PURITY);
        MinSomaticPuritySpread = configBuilder.getDecimal(SOMATIC_MIN_PURITY_SPREAD);
        SomaticPenaltyWeight = configBuilder.getDecimal(SOMATIC_PENALTY_WEIGHT);
        HighlyDiploidPercentage = configBuilder.getDecimal(HIGHLY_DIPLOID_PERCENTAGE);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(
                SOMATIC_MIN_PEAK, "Minimum number of somatic variants to consider a peak", SOMATIC_MIN_PEAK_DEFAULT);

        configBuilder.addInteger(SOMATIC_MIN_TOTAL,
                "Minimum number of somatic variants required to assist highly diploid fits", SOMATIC_MIN_VARIANTS_DEFAULT);

        configBuilder.addDecimal(
                SOMATIC_MIN_PURITY, "Minimum spread within candidate purities before somatics can be used",
                SOMATIC_MIN_PURITY_DEFAULT);

        configBuilder.addDecimal(
                SOMATIC_MIN_PURITY_SPREAD, "Minimum spread within candidate purities before somatics can be used",
                SOMATIC_MIN_PURITY_SPREAD_DEFAULT);

        configBuilder.addDecimal(
                SOMATIC_PENALTY_WEIGHT, "Proportion of somatic deviation to include in fitted purity score",
                SOMATIC_PENALTY_WEIGHT_DEFAULT);

        configBuilder.addDecimal(
                HIGHLY_DIPLOID_PERCENTAGE, "Proportion of genome that must be diploid before using somatic fit",
                HIGHLY_DIPLOID_PERCENTAGE_DEFAULT);
    }
}
