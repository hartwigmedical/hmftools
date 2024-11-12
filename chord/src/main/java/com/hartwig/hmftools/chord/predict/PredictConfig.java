package com.hartwig.hmftools.chord.predict;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PredictConfig
{
    public final String MutContextsFile;
    public final String OutputFile;

    private static final String MUT_CONTEXTS_FILE = "mut_contexts_file";
    private static final String MUT_CONTEXTS_FILE_DESC = "Path to *.chord.mutation_contexts.tsv file";

    private static final String OUTPUT_FILE = "output_file";
    private static final String OUTPUT_FILE_DESC = "Path to output predictions";

    public PredictConfig(final ConfigBuilder configBuilder)
    {
        MutContextsFile = configBuilder.getValue(MUT_CONTEXTS_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(MUT_CONTEXTS_FILE, true, MUT_CONTEXTS_FILE_DESC);
        configBuilder.addConfigItem(OUTPUT_FILE, true, OUTPUT_FILE_DESC);
    }
}
