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
import static com.hartwig.hmftools.neo.scorer.SampleData.loadFromConfig;

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
    public final String CohortSampleTpmFile;
    public final String CohortTpmMediansFile;

    public final List<SampleData> Samples;

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
    public static final String COHORT_SAMPLE_TPM_FILE = "cohort_trans_exp_file";
    public static final String COHORT_TPM_MEDIANS_FILE = "cancer_tpm_medians_file";

    private static final String LIKELIHOOD_THRESHOLD = "rank_threshold";
    private static final String SIMILARITY_THRESHOLD = "sim_threshold";
    public static final String PEPTIDE_LENGTHS = "peptide_lengths";
    public static final String PEPTIDE_FLANKS = "peptide_flanks";

    private static final String WRITE_TYPES = "write_types";

    public static final String RNA_SAMPLE_APPEND_SUFFIX = "_RNA";

    public NeoScorerConfig(final CommandLine cmd)
    {
        Samples = loadFromConfig(cmd);

        String sampleDataDir = cmd.hasOption(SAMPLE_DATA_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR)) : "";

        NeoDataDir = cmd.getOptionValue(NEO_DATA_DIR, sampleDataDir);
        LilacDataDir = cmd.getOptionValue(LILAC_DATA_DIR, sampleDataDir);
        IsofoxDataDir = cmd.getOptionValue(ISF_DATA_DIR, sampleDataDir);
        RnaSomaticVcf = cmd.getOptionValue(RNA_SOMATIC_VCF, sampleDataDir);
        OutputDir = cmd.hasOption(OUTPUT_DIR) ? parseOutputDir(cmd) : sampleDataDir;
        RnaSampleSuffix = cmd.getOptionValue(RNA_SAMPLE_SUFFIX, RNA_SAMPLE_APPEND_SUFFIX);

        CohortSampleTpmFile = cmd.getOptionValue(COHORT_SAMPLE_TPM_FILE);
        CohortTpmMediansFile = cmd.getOptionValue(COHORT_TPM_MEDIANS_FILE);

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

        if(Samples.size() == 1)
            filename += Samples.get(0).Id + ".neo";
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
        options.addOption(COHORT_SAMPLE_TPM_FILE, true, "Cohort gene expression matrix");
        options.addOption(COHORT_TPM_MEDIANS_FILE, true, "TPM medians per cancer type and pan-cancer");

        options.addOption(PEPTIDE_LENGTHS, true, "Peptide length min-max, separated by '-', eg 8-12");
        options.addOption(PEPTIDE_FLANKS, true, "Peptide flanking amino acids");

        ScoreConfig.addCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        options.addOption(LIKELIHOOD_THRESHOLD, true, "Rank threshold to write full peptide data");
        options.addOption(SIMILARITY_THRESHOLD, true, "Immunogenic similarity threshold to write full peptide data");

        options.addOption(WRITE_TYPES, true, "Valid types: ALLELE_PEPTIDE, PEPTIDE, NEOEPITOPE, SAMPLE_SUMMARY sep by ';'");
        addOutputOptions(options);
        addThreadOptions(options);
    }
}
