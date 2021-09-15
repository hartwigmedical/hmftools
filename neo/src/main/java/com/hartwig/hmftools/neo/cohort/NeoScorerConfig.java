package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.NeoCommon.THREADS;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleIdsFile;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class NeoScorerConfig
{
    public final String NeoDataDir;
    public final String LilacDataDir;
    public final String IsofoxDataDir;
    public final String SampleTranscriptExpressionFile;

    public final List<String> SampleIds;

    public final String OutputDir;
    public final String OutputId;
    public final List<OutputType> WriteTypes;

    public final double LikelihoodThreshold;
    public final int Threads;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String LILAC_DATA_DIR = "lilac_data_dir";
    public static final String PREDICTION_DATA_DIR = "mcf_prediction_dir";
    public static final String ISF_DATA_DIR = "isofox_neo_dir";
    public static final String SAMPLE_TRANS_EXP_FILE = "sample_trans_exp_file";

    private static final String LIKELIHOOD_THRESHOLD = "rank_threshold";

    private static final String WRITE_TYPES = "write_types";

    public NeoScorerConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), SampleIds);

        NeoDataDir = cmd.getOptionValue(NEO_DATA_DIR);
        LilacDataDir = cmd.getOptionValue(LILAC_DATA_DIR);
        IsofoxDataDir = cmd.getOptionValue(ISF_DATA_DIR);

        SampleTranscriptExpressionFile = cmd.getOptionValue(SAMPLE_TRANS_EXP_FILE);

        LikelihoodThreshold = Double.parseDouble(cmd.getOptionValue(LIKELIHOOD_THRESHOLD, "0.02"));

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        WriteTypes = Lists.newArrayList();

        if(cmd.hasOption(WRITE_TYPES))
        {
            String[] types = cmd.getOptionValue(WRITE_TYPES).split(";");
            Arrays.stream(types).map(x -> OutputType.valueOf(x)).forEach(x -> WriteTypes.add(x));
        }
    }

    public String formFilename(final String fileId)
    {
        if(OutputId == null || OutputId.isEmpty())
            return String.format("%sneo_cohort_%s.csv", OutputDir, fileId);
        else
            return String.format("%sneo_cohort_%s_%s.csv", OutputDir, OutputId, fileId);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(LILAC_DATA_DIR, true, "Directory for Lilac coverage files");
        options.addOption(ISF_DATA_DIR, true, "Directory for Isofox neoepitope coverage files");
        options.addOption(SAMPLE_TRANS_EXP_FILE, true, "Cohort gene expression matrix");

        ScoreConfig.addCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        options.addOption(LIKELIHOOD_THRESHOLD, true, "Rank threshold to write full peptide data");
        options.addOption(THREADS, true, "Thread count");

        options.addOption(WRITE_TYPES, true, "Valid types: ALLELE_PEPTIDE, PEPTIDE, NEOEPITOPE, SAMPLE_SUMMARY sep by ';'");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file ID");
    }
}
