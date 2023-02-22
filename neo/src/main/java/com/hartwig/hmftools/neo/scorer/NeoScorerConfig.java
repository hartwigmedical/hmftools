package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

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
    public final String RnaSomaticVcf;
    public final String IsofoxDataDir;
    public final String CohortTranscriptExpressionFile;

    public final List<String> SampleIds;

    public final String OutputDir;
    public final String OutputId;
    public final String RnaSampleSuffix;
    public final List<OutputType> WriteTypes;

    public final double LikelihoodThreshold;
    public final double SimilarityThreshold;
    public final int Threads;

    public static final String SAMPLE = "sample";
    public static final String SAMPLE_DATA_DIR = "sample_data_dir";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String LILAC_DATA_DIR = "lilac_data_dir";
    public static final String RNA_SOMATIC_VCF = "rna_somatic_vcf";
    public static final String RNA_SAMPLE_SUFFIX = "rna_sample_suffix";
    public static final String PREDICTION_DATA_DIR = "mcf_prediction_dir";
    public static final String ISF_DATA_DIR = "isofox_neo_dir";
    public static final String COHORT_TRANS_EXP_FILE = "cohort_trans_exp_file";

    private static final String LIKELIHOOD_THRESHOLD = "rank_threshold";
    private static final String SIMILARITY_THRESHOLD = "sim_threshold";

    private static final String WRITE_TYPES = "write_types";

    public static final String RNA_SAMPLE_APPEND_SUFFIX = "_RNA";

    public NeoScorerConfig(final CommandLine cmd)
    {
        if(cmd.hasOption(SAMPLE))
            SampleIds = Lists.newArrayList(cmd.getOptionValue(SAMPLE));
        else
            SampleIds = loadSampleIdsFile(cmd);

        String sampleDataDir = cmd.hasOption(SAMPLE_DATA_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR)) : "";

        NeoDataDir = cmd.getOptionValue(NEO_DATA_DIR, sampleDataDir);
        LilacDataDir = cmd.getOptionValue(LILAC_DATA_DIR, sampleDataDir);
        IsofoxDataDir = cmd.getOptionValue(ISF_DATA_DIR, sampleDataDir);
        RnaSomaticVcf = cmd.getOptionValue(RNA_SOMATIC_VCF, sampleDataDir);
        OutputDir = cmd.hasOption(OUTPUT_DIR) ? parseOutputDir(cmd) : sampleDataDir;
        RnaSampleSuffix = cmd.getOptionValue(RNA_SAMPLE_SUFFIX, RNA_SAMPLE_APPEND_SUFFIX);

        CohortTranscriptExpressionFile = cmd.getOptionValue(COHORT_TRANS_EXP_FILE);

        LikelihoodThreshold = Double.parseDouble(cmd.getOptionValue(LIKELIHOOD_THRESHOLD, "0.02"));
        SimilarityThreshold = Double.parseDouble(cmd.getOptionValue(SIMILARITY_THRESHOLD, "10"));

        OutputId = cmd.getOptionValue(OUTPUT_ID);
        Threads = parseThreads(cmd);

        WriteTypes = Lists.newArrayList();

        if(cmd.hasOption(WRITE_TYPES))
        {
            String[] types = cmd.getOptionValue(WRITE_TYPES).split(";");
            Arrays.stream(types).map(x -> OutputType.valueOf(x)).forEach(x -> WriteTypes.add(x));
        }
    }

    public String formFilename(final String fileId)
    {
        String filename = OutputDir;

        if(SampleIds.size() == 1)
            filename += SampleIds.get(0) + ".neo";
        else
            filename += "neo_cohort";

        if(OutputId != null && !OutputId.isEmpty())
            filename += "." + OutputId;

        filename += "." + fileId;
        filename += ".csv";

        return filename;
    }

    public static void addCmdLineArgs(Options options)
    {
        addSampleIdFile(options);
        options.addOption(SAMPLE, true, "Sample ID for single sample");
        options.addOption(SAMPLE_DATA_DIR, true, "Directory for sample files");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(LILAC_DATA_DIR, true, "Directory for Lilac coverage files");
        options.addOption(RNA_SOMATIC_VCF, true, "Directory for Purple somatic variant RNA-appended files");
        options.addOption(RNA_SAMPLE_SUFFIX, true, "RNA sample suffix in Sage-appended VCF");
        options.addOption(ISF_DATA_DIR, true, "Directory for Isofox neoepitope coverage files");
        options.addOption(COHORT_TRANS_EXP_FILE, true, "Cohort gene expression matrix");

        ScoreConfig.addCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        options.addOption(LIKELIHOOD_THRESHOLD, true, "Rank threshold to write full peptide data");
        options.addOption(SIMILARITY_THRESHOLD, true, "Immunogenic similarity threshold to write full peptide data");

        options.addOption(WRITE_TYPES, true, "Valid types: ALLELE_PEPTIDE, PEPTIDE, NEOEPITOPE, SAMPLE_SUMMARY sep by ';'");
        addOutputOptions(options);
        addThreadOptions(options);
    }
}
