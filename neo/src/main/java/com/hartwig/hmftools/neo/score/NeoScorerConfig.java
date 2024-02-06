package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.NEO_FILE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.NEO_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.NEO_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.neo.score.SampleData.loadFromConfig;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

public class NeoScorerConfig
{
    public final String NeoDir;
    public final String LilacDir;
    public final String RnaSomaticVcf;
    public final String IsofoxDir;
    public final String PurpleDir;
    public final String CohortSampleTpmFile;
    public final String CohortTpmMediansFile;

    public final List<SampleData> Samples;

    public final String OutputDir;
    public final String OutputId;
    public final String RnaSampleId;
    public final List<OutputType> WriteTypes;

    public final double LikelihoodThreshold;
    public final double SimilarityThreshold;
    public final int Threads;

    public static final String CANCER_TYPE = "cancer_type";
    public static final String RNA_SOMATIC_VCF = "rna_somatic_vcf";
    public static final String COHORT_SAMPLE_TPM_FILE = "cohort_trans_exp_file";

    public static final String COHORT_TPM_MEDIANS_FILE = "cancer_tpm_medians_file";
    public static final String COHORT_TPM_MEDIANS_FILE_DESC = "TPM medians per cancer type and pan-cancer";

    public static final String LIKELIHOOD_THRESHOLD = "rank_threshold";
    private static final String SIMILARITY_THRESHOLD = "sim_threshold";
    public static final String PEPTIDE_LENGTHS = "peptide_lengths";

    private static final String WRITE_TYPES = "write_types";

    public static final String RNA_SAMPLE_ID_SUFFIX = "_RNA";

    public NeoScorerConfig(final ConfigBuilder configBuilder)
    {
        Samples = loadFromConfig(configBuilder);

        String sampleDataDir = configBuilder.hasValue(SAMPLE_DATA_DIR_CFG) ? checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG)) : "";

        NeoDir = checkAddDirSeparator(configBuilder.getValue(NEO_DIR_CFG, sampleDataDir));
        LilacDir = checkAddDirSeparator(configBuilder.getValue(LILAC_DIR_CFG, sampleDataDir));
        PurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG, sampleDataDir));
        IsofoxDir = checkAddDirSeparator(configBuilder.getValue(ISOFOX_DIR_CFG, sampleDataDir));
        RnaSomaticVcf = configBuilder.getValue(RNA_SOMATIC_VCF, sampleDataDir);
        OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : sampleDataDir;
        RnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);

        CohortSampleTpmFile = configBuilder.getValue(COHORT_SAMPLE_TPM_FILE);
        CohortTpmMediansFile = configBuilder.getValue(COHORT_TPM_MEDIANS_FILE);

        LikelihoodThreshold = configBuilder.getDecimal(LIKELIHOOD_THRESHOLD);
        SimilarityThreshold = configBuilder.getDecimal(SIMILARITY_THRESHOLD);

        OutputId = configBuilder.getValue(OUTPUT_ID);
        Threads = parseThreads(configBuilder);

        WriteTypes = Lists.newArrayList();

        if(configBuilder.hasValue(WRITE_TYPES))
        {
            String[] types = configBuilder.getValue(WRITE_TYPES).split(ITEM_DELIM);
            Arrays.stream(types).map(x -> OutputType.valueOf(x)).forEach(x -> WriteTypes.add(x));
        }
        else
        {
            WriteTypes.add(OutputType.NEOEPITOPE);
            WriteTypes.add(OutputType.ALLELE_PEPTIDE);
        }
    }

    public String formFilename(final String fileId)
    {
        String filename = OutputDir;

        if(Samples.size() == 1)
            filename += Samples.get(0).TumorId + "." + NEO_FILE_ID;
        else
            filename += "neo_cohort";

        if(OutputId != null && !OutputId.isEmpty())
            filename += "." + OutputId;

        filename += "." + fileId;
        filename += TSV_EXTENSION;

        return filename;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addConfigItem(CANCER_TYPE, false, "Cancer type for sample, matching those in cohort TPM medians file");
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);
        configBuilder.addPath(NEO_DIR_CFG, false, NEO_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(LILAC_DIR_CFG, false, LILAC_DIR_DESC);
        configBuilder.addPath(RNA_SOMATIC_VCF, false, "Directory for Purple somatic variant RNA-appended files");
        configBuilder.addConfigItem(RNA_SAMPLE_ID, RNA_SAMPLE_DESC + " in Sage-appended VCF");
        configBuilder.addPath(ISOFOX_DIR_CFG, false, "Directory for Isofox files (Transcript expresion, Neoepitope coverage)");
        configBuilder.addPath(COHORT_SAMPLE_TPM_FILE, false, "Cohort gene expression matrix");
        configBuilder.addPath(COHORT_TPM_MEDIANS_FILE, false, COHORT_TPM_MEDIANS_FILE_DESC);

        configBuilder.addConfigItem(PEPTIDE_LENGTHS, false, "Peptide length min-max, separated by '-', eg 8-12");

        ScoreConfig.registerConfig(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        EnsemblDataCache.addEnsemblDir(configBuilder);

        configBuilder.addDecimal(
                LIKELIHOOD_THRESHOLD, "Rank threshold to write full peptide data, default 0 (not applied)", 0);

        configBuilder.addDecimal(
                SIMILARITY_THRESHOLD, "Immunogenic similarity threshold to write full peptide data", 10);

        configBuilder.addConfigItem(
                WRITE_TYPES, false, "Valid types: ALLELE_PEPTIDE, PEPTIDE, NEOEPITOPE, SAMPLE_SUMMARY sep by ';'");

        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
