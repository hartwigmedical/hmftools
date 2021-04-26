package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class SomaticFitConfig
{
    public static final String SOMATIC_MIN_PEAK = "somatic_min_peak";
    public static final String SOMATIC_MIN_TOTAL = "somatic_min_variants";
    public static final String SOMATIC_MIN_PURITY = "somatic_min_purity";
    public static final String SOMATIC_MIN_PURITY_SPREAD = "somatic_min_purity_spread";
    public static final String SOMATIC_PENALTY_WEIGHT = "somatic_penalty_weight";
    public static final String HIGHLY_DIPLOID_PERCENTAGE = "highly_diploid_percentage";

    public static final double SOMATIC_MIN_PURITY_DEFAULT = 0.17;
    public static final double SOMATIC_MIN_PURITY_SPREAD_DEFAULT = 0.15;
    public static final int SOMATIC_MIN_PEAK_DEFAULT = 10;
    public static final int SOMATIC_MIN_VARIANTS_DEFAULT = 10;
    public static final double SOMATIC_PENALTY_WEIGHT_DEFAULT = 1;
    public static final double HIGHLY_DIPLOID_PERCENTAGE_DEFAULT = 0.97;

    public int MinTotalVariants;
    public int MinPeakVariants;
    public double MinSomaticPurity;
    public double MinSomaticPuritySpread;
    public double SomaticPenaltyWeight;
    public double HighlyDiploidPercentage;

    public int minTotalSomaticVariantAlleleReadCount()
    {
        return 5000;
    }
    public int minTotalSvFragmentCount()
    {
        return 1000;
    }
    public double clonalityMaxPloidy()
    {
        return 10;
    }
    public double clonalityBinWidth()
    {
        return 0.05;
    }

    public SomaticFitConfig(@NotNull CommandLine cmd)
    {
        MinTotalVariants = getConfigValue(cmd, SOMATIC_MIN_TOTAL, SOMATIC_MIN_VARIANTS_DEFAULT);
        MinPeakVariants = getConfigValue(cmd, SOMATIC_MIN_PEAK, SOMATIC_MIN_PEAK_DEFAULT);
        MinSomaticPurity = getConfigValue(cmd, SOMATIC_MIN_PURITY, SOMATIC_MIN_PURITY_DEFAULT);
        MinSomaticPuritySpread = getConfigValue(cmd, SOMATIC_MIN_PURITY_SPREAD, SOMATIC_MIN_PURITY_SPREAD_DEFAULT);
        SomaticPenaltyWeight = getConfigValue(cmd, SOMATIC_PENALTY_WEIGHT, SOMATIC_PENALTY_WEIGHT_DEFAULT);
        HighlyDiploidPercentage = getConfigValue(cmd, HIGHLY_DIPLOID_PERCENTAGE, HIGHLY_DIPLOID_PERCENTAGE_DEFAULT);
    }


    public static void addOptions(@NotNull Options options)
    {
        options.addOption(SOMATIC_MIN_PEAK, true,
                "Minimum number of somatic variants to consider a peak. Default " + SOMATIC_MIN_PEAK_DEFAULT);
        options.addOption(SOMATIC_MIN_TOTAL, true,
                "Minimum number of somatic variants required to assist highly diploid fits. Default " + SOMATIC_MIN_VARIANTS_DEFAULT);
        options.addOption(SOMATIC_MIN_PURITY, true,
                "Somatic fit will not be used if both somatic and fitted purities are less than this value. Default "
                        + SOMATIC_MIN_PURITY_DEFAULT);
        options.addOption(SOMATIC_MIN_PURITY_SPREAD, true,
                "Minimum spread within candidate purities before somatics can be used. Default " + SOMATIC_MIN_PURITY_SPREAD_DEFAULT);
        options.addOption(SOMATIC_PENALTY_WEIGHT, true,
                "Proportion of somatic deviation to include in fitted purity score. Default " + SOMATIC_PENALTY_WEIGHT_DEFAULT);
        options.addOption(HIGHLY_DIPLOID_PERCENTAGE, true,
                "Proportion of genome that must be diploid before using somatic fit. Default " + HIGHLY_DIPLOID_PERCENTAGE_DEFAULT);
    }
}
