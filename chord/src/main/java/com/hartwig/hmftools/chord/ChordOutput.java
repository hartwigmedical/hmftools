package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;

import com.hartwig.hmftools.chord.ChordConfig;

public class ChordOutput
{
    public static final String COHORT_FILE_PREFIX = "cohort";

    public static final String CHORD_MUTATION_CONTEXTS = ".chord.mutation_contexts";
    public static final String CHORD_PREDICTIONS = ".chord.predictions";

    private static String formOutputFile(String outputDir, String prefix, String suffix, String outputId)
    {
        String outputFile = outputDir + "/" + prefix;

        outputFile += suffix;

        if(outputId != null && !outputId.isEmpty())
            outputFile += "." + outputId;

        outputFile += TSV_EXTENSION;

        return outputFile;
    }

    private static String prefixSampleIdOrCohort(ChordConfig config)
    {
        return config.isSingleSample() ?
                config.SampleIds.get(0) :
                COHORT_FILE_PREFIX;
    }

    public static String mutationContextsFile(ChordConfig config)
    {
        return formOutputFile(config.OutputDir, prefixSampleIdOrCohort(config), CHORD_MUTATION_CONTEXTS, config.OutputId);
    }

    public static String predictionsFile(ChordConfig config)
    {
        return formOutputFile(config.OutputDir, prefixSampleIdOrCohort(config), CHORD_PREDICTIONS, config.OutputId);
    }
}
