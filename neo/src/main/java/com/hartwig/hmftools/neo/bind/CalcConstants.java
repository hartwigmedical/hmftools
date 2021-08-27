package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT_MAX;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_FACTOR;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEP_LEN_WEIGHT_MAX;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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

    public CalcConstants(final CommandLine cmd)
    {
        PeptideLengthWeightFactor = parseParam(cmd, PEP_LEN_WEIGHT_FACTOR, DEFAULT_PEP_LEN_WEIGHT_FACTOR);

        PeptideLengthWeightMax = parseParam(cmd, PEP_LEN_WEIGHT_MAX, DEFAULT_PEP_LEN_WEIGHT_MAX);

        AlleleWeightFactor = parseParam(cmd, ALLELE_MOTIF_WEIGHT_FACTOR, DEFAULT_ALLELE_MOTIF_WEIGHT_FACTOR);

        AlleleWeightMax = parseParam(cmd, ALLELE_MOTIF_WEIGHT_MAX, DEFAULT_ALLELE_MOTIF_WEIGHT_MAX);
        NoiseProbability = parseParam(cmd, NOISE_PROB, 0); // DEFAULT_NOISE_PROB)
        NoiseWeight = parseParam(cmd, NOISE_WEIGHT, 0); // DEFAULT_NOISE_WEIGHT)
    }

    private static double parseParam(final CommandLine cmd, final String config, double defaultValue)
    {
        return Double.parseDouble(cmd.getOptionValue(config, String.valueOf(defaultValue)));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(ALLELE_MOTIF_WEIGHT_MAX, true, "Allele motif weight");
        options.addOption(PEP_LEN_WEIGHT_FACTOR, true, "Length weight");
        options.addOption(NOISE_PROB, true, "Noise target probability");
        options.addOption(NOISE_WEIGHT, true, "Noise weight");
    }

}
