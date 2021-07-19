package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.immutables.value.internal.$guava$.collect.$MutableClassToInstanceMap;

public class BinderConfig
{
    public final String TrainingDataFile;
    public final String OutputDir;
    public final String OutputId;

    public final List<String> SpecificAlleles;

    public final List<double[]> BindingLevelScores;

    private static final String TRAINING_DATA_FILE = "training_data_file";
    private static final String SPECIFIC_ALLELES = "specific_alleles";
    private static final String OUTPUT_ID = "output_id";

    private static final String BIND_LEVEL_SCORES = "binding_level_scores";

    public BinderConfig(final CommandLine cmd)
    {
        TrainingDataFile = cmd.getOptionValue(TRAINING_DATA_FILE);
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        BindingLevelScores = Lists.newArrayList();

        String[] levelScores = cmd.getOptionValue(BIND_LEVEL_SCORES).split(ITEM_DELIM, -1);

        for(String levelScore :levelScores)
        {
            String[] item = levelScore.split("=", -1);
            int level = Integer.parseInt(item[0]);
            double score = Double.parseDouble(item[1]);
            BindingLevelScores.add(new double[] {level, score});
        }

        SpecificAlleles = Lists.newArrayList();

        if(cmd.hasOption(SPECIFIC_ALLELES))
        {
            Arrays.stream(cmd.getOptionValue(SPECIFIC_ALLELES).split(ITEM_DELIM, -1)).forEach(x -> SpecificAlleles.add(x));
            NE_LOGGER.info("filtering for {} alleles: {}", SpecificAlleles.size(), SpecificAlleles);
        }
    }

    public String formFilename(final String fileId)
    {
        if(OutputId.isEmpty())
            return String.format("%sbind_%s.csv", OutputDir, fileId);
        else
            return String.format("%sbind_%s_%s.csv", OutputDir, OutputId, fileId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(TRAINING_DATA_FILE, true, "Directory for sample prediction result files");
        options.addOption(SPECIFIC_ALLELES, true, "List of alleles separated by ';'");
        options.addOption(BIND_LEVEL_SCORES, true, "Affinities below the level deemed to bind");

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }

}
