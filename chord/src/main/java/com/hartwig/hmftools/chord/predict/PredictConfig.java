package com.hartwig.hmftools.chord.predict;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PredictConfig
{
    public final String MutContextsFile;
    public final String OutputFile;
    public final String ChordModelFile;

    private static final String MUT_CONTEXTS_FILE = "mut_contexts_file";
    private static final String MUT_CONTEXTS_FILE_DESC = "Path to *.chord.mutation_contexts.tsv file";

    private static final String OUTPUT_FILE = "output_file";
    private static final String OUTPUT_FILE_DESC = "Path to output predictions";

    private static final String CHORD_MODEL_FILE = "chord_model_file";
    private static final String CHORD_MODEL_FILE_DESC = "Path to the CHORD.rds file";

    public PredictConfig(final ConfigBuilder configBuilder)
    {
        MutContextsFile = configBuilder.getValue(MUT_CONTEXTS_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        ChordModelFile = configBuilder.getValue(CHORD_MODEL_FILE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(MUT_CONTEXTS_FILE, true, MUT_CONTEXTS_FILE_DESC);
        configBuilder.addConfigItem(OUTPUT_FILE, true, OUTPUT_FILE_DESC);
        configBuilder.addPath(CHORD_MODEL_FILE, false, CHORD_MODEL_FILE_DESC);
    }
}
