package com.hartwig.hmftools.data_analyser.calcs;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SimConfig {

    final public int SampleCount;
    final public String SigFactorsFilename;
    final public String SignaturesFilename;
    final public boolean SeedRandom;
    final public boolean ApplyNoise;

    public static String SIG_FACTOR_FILENAME = "sig_factors_file";
    public static String SIGNATURES_FILENAME = "signature_file";
    public static String SIM_SAMPLE_COUNT = "sim_sample_count";
    public static String SIM_RANDOM_SEED = "sim_seed_random";
    public static String SIM_APPLY_NOISE = "sim_apply_noise";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SIGNATURES_FILENAME, true, "File containing signature defintions");
        options.addOption(SIG_FACTOR_FILENAME, true, "File containing signature parameters");
        options.addOption(SIM_SAMPLE_COUNT, true, "Number of samples in cohort");
        options.addOption(SIM_RANDOM_SEED, false, "Whether to seed random number generation");
        options.addOption(SIM_APPLY_NOISE, false, "Apply poisson noise to bucket counts");
    }

    public SimConfig(final CommandLine cmd)
    {
        SampleCount = Integer.parseInt(cmd.getOptionValue(SIM_SAMPLE_COUNT));
        SignaturesFilename = cmd.getOptionValue(SIGNATURES_FILENAME);
        SigFactorsFilename = cmd.getOptionValue(SIG_FACTOR_FILENAME);
        SeedRandom = cmd.hasOption(SIM_RANDOM_SEED);
        ApplyNoise = cmd.hasOption(SIM_APPLY_NOISE);
    }
}
