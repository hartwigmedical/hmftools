package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_MAX;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_MAX;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class CalcConstants
{
    public final double PeptideLengthWeightFactor;
    public final double PeptideLengthWeightMax;

    public final double AlleleWeightFactor;
    public final double AlleleWeightMax;

    public final double NoiseProbability;
    public final double NoiseWeight;

    private static final String PEP_LEN_WEIGHT_FACTOR = "pl_weight_factor";
    private static final String PEP_LEN_WEIGHT_MAX = "pl_weight_max";
    private static final String ALLELE_MOTIF_WEIGHT_FACTOR = "am_weight_factor";
    private static final String ALLELE_MOTIF_WEIGHT_MAX = "am_weight_max";
    private static final String NOISE_PROB = "noise_prob";
    private static final String NOISE_WEIGHT = "noise_weight";

    public CalcConstants(
            double peptideLengthWeightFactor, double peptideLengthWeightMax, double alleleWeightFactor, double alleleWeightMax,
            double noiseProb, double noiseWeight)
    {
        PeptideLengthWeightFactor = peptideLengthWeightFactor;
        PeptideLengthWeightMax = peptideLengthWeightMax;
        AlleleWeightFactor = alleleWeightFactor;
        AlleleWeightMax = alleleWeightMax;
        NoiseProbability = noiseProb;
        NoiseWeight = noiseWeight;
    }

    public CalcConstants(final ConfigBuilder configBuilder)
    {
        PeptideLengthWeightFactor = configBuilder.getDecimal(PEP_LEN_WEIGHT_FACTOR);
        PeptideLengthWeightMax = configBuilder.getDecimal(PEP_LEN_WEIGHT_MAX);
        AlleleWeightFactor = configBuilder.getDecimal(ALLELE_MOTIF_WEIGHT_FACTOR);
        AlleleWeightMax = configBuilder.getDecimal(ALLELE_MOTIF_WEIGHT_MAX);
        NoiseProbability = configBuilder.getDecimal(NOISE_PROB);
        NoiseWeight = configBuilder.getDecimal(NOISE_WEIGHT);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimalItem(PEP_LEN_WEIGHT_FACTOR, false, "Length weight factor", DEFAULT_PEP_LEN_WEIGHT_FACTOR);
        configBuilder.addDecimalItem(PEP_LEN_WEIGHT_MAX, false, "Length weight max", DEFAULT_PEP_LEN_WEIGHT_MAX);

        configBuilder.addDecimalItem(
                ALLELE_MOTIF_WEIGHT_FACTOR, false, "Allele motif weight factor", DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR);

        configBuilder.addDecimalItem(
                ALLELE_MOTIF_WEIGHT_MAX, false, "Allele motif weight max", DEFAULT_ALLELE_MOTIF_WEIGHT_MAX);

        configBuilder.addDecimalItem(NOISE_PROB, false, "Noise target probability", 0);

        configBuilder.addDecimalItem(NOISE_WEIGHT, false, "Noise weight", 0);
    }

}
