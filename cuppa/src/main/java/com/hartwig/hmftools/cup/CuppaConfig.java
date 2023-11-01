package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.cup.common.CupConstants.COMBINED_DAMPEN_FACTOR_DEFAULT;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FEATURE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_SIG_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_SV_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_TRAITS_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_CANCER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_CANCER_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_DRIVER_AVG;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_FEATURE_PREV;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENDER_RATES;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_CANCER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SIG_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_SIGNATURES;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SV_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_PERC;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_TRAIT_RATES;
import static com.hartwig.hmftools.common.cuppa.CategoryType.COMBINED;
import static com.hartwig.hmftools.common.cuppa.CategoryType.isDna;
import static com.hartwig.hmftools.common.cuppa.CategoryType.isRna;
import static com.hartwig.hmftools.cup.common.CupConstants.DEFAULT_RNA_LENGTH;
import static com.hartwig.hmftools.cup.common.CupConstants.FEATURE_DAMPEN_FACTOR_DEFAULT;
import static com.hartwig.hmftools.cup.prep.PrepConfig.addPipelineDirectories;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.REF_COHORT_FEATURES_FILE;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.REF_COHORT_SAMPLE_TRAITS_FILE;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.REF_COHORT_SIG_CONTRIBS_FILE;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.REF_COHORT_SV_DATA_FILE;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.rna.AltSjClassifier;
import com.hartwig.hmftools.cup.rna.GeneExpressionClassifier;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;
import com.hartwig.hmftools.cup.traits.SampleTraitClassifier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CuppaConfig
{
    public final RefGenomeVersion RefGenVersion;
    public final List<CategoryType> Categories;

    // reference data
    public final String RefDataDir;

    public final String RefSampleDataFile;
    public final String RefSnvCountsFile;
    public final String RefSvPercFile;
    public final String RefSigContributionFile;
    public final String RefFeaturePrevFile;
    public final String RefTraitPercFile;
    public final String RefTraitRateFile;
    public final String RefGenderRateFile;
    public final String RefSnvCancerPosFreqFile;
    public final String RefSnvSamplePosFreqFile;
    public final String RefDriverAvgFile;

    public final String RefGeneExpCancerFile;
    public final String RefGeneExpSampleFile;
    public final String RefAltSjCancerFile;
    public final String RefAltSjSampleFile;
    public final String RefSnvSignaturesFile;

    // a single sample directory
    public final String SampleDataDir;

    // or standard pipeline directories
    public final String LinxDir;
    public final String PurpleDir;
    public final String VirusDir;
    public final String IsofoxDir;

    // cohort files, formed during ref data building
    public final String SampleDataFile;
    public final String RefSampleFeatureFile;
    public final String RefSampleTraitsFile;
    public final String RefSampleSigContribFile;
    public final String RefSampleSvFile;

    public final boolean TestRefData;

    public final NoiseRefCache NoiseAdjustments;
    public final boolean NoSubtypeCollapse;
    public final double FeatureDampenFactor;
    public final double CombinedDampenFactor;

    // database access
    public final DatabaseAccess DbAccess;

    public final boolean WriteSimilarities;
    public final boolean WriteDetailedScores;
    public final boolean WriteCondensed;
    public final boolean CreatePdf;

    public final String OutputDir;
    public final String OutputFileId;
    public final int Threads;

    // config strings
    public static final String CATEGORIES = "categories";

    public static final String ALL_CATEGORIES = "ALL";
    public static final String DNA_CATEGORIES = "DNA";
    public static final String RNA_CATEGORIES = "RNA";

    // either a single sample to be tested for a file containing the samples to be tested
    public static final String SAMPLE_RNA_LENGTH = "sample_rna_length";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";

    public static final String REF_DATA_DIR = "ref_data_dir";
    public static final String TEST_REF_SAMPLE_DATA = "test_ref_sample_data";

    // reference data files
    public static final String REF_SAMPLE_DATA_FILE = "ref_sample_data_file";
    public static final String REF_SNV_COUNTS_FILE = "ref_snv_counts_file";
    public static final String REF_SNV_SAMPLE_POS_FREQ_FILE = "ref_sample_snv_pos_freq_file";
    private static final String REF_SNV_CANCER_POS_FREQ_FILE = "ref_cancer_snv_pos_freq_file";
    private static final String REF_SIG_CONTRIB_FILE = "ref_sig_contrib_file";
    private static final String REF_FEAT_PREV_FILE = "ref_feature_prev_file";
    private static final String REF_DRIVER_AVG_FILE = "ref_feature_avg_file";
    private static final String REF_TRAIT_PERC_FILE = "ref_trait_perc_file";
    private static final String REF_TRAIT_RATE_FILE = "ref_trait_rate_file";
    private static final String REF_GENDER_RATE_FILE = "ref_gender_rate_file";
    private static final String REF_SV_PERC_FILE = "ref_sv_perc_file";
    private static final String REF_RNA_GENE_EXP_CANCER_FILE = "ref_gene_exp_cancer_file";
    public static final String REF_RNA_GENE_EXP_SAMPLE_FILE = "ref_gene_exp_sample_file";
    private static final String REF_RNA_ALT_SJ_CANCER_FILE = "ref_alt_sj_cancer_file";
    public static final String REF_RNA_ALT_SJ_SAMPLE_FILE = "ref_alt_sj_sample_file";
    public static final String REF_SNV_SIGNATURES_FILE = "ref_snv_signatures_file";

    public static final String NOISE_ALLOCATIONS = "noise_allocations";
    public static final String NOISE_ALLOCATIONS_DESC = "Noise allocations by classifier type, or 'NONE' or 'DEFAULTS'";
    public static final String NO_SUBTYPE_COLLAPSE = "no_subtype_collapse";

    public static final String FEATURE_DAMPEN_FACTOR = "feature_dampen_factor";
    public static final String COMBINED_DAMPEN_FACTOR = "combined_dampen_factor";

    public static final String WRITE_SIMS = "write_similarities";
    public static final String WRITE_DETAILED_SCORES = "write_detailed_scores";
    public static final String WRITE_CONDENSED = "write_condensed";
    public static final String CREATE_PDF = "create_pdf";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CuppaConfig.class);

    // file fields
    public static final String FLD_SAMPLE_ID = "SampleId";
    public static final String FLD_CANCER_TYPE = "CancerType";
    public static final String FLD_CANCER_SUBTYPE = "CancerSubtype";
    public static final String FLD_RNA_READ_LENGTH = "RnaReadLength";
    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";

    public CuppaConfig(final ConfigBuilder configBuilder)
    {
        Categories = configCategories(configBuilder);

        CUP_LOGGER.info("running classifiers: {}", Categories.isEmpty() ? ALL_CATEGORIES : Categories.toString());

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        RefDataDir = checkAddDirSeparator(configBuilder.getValue(REF_DATA_DIR, ""));

        RefSampleDataFile = getRefDataFile(configBuilder, REF_SAMPLE_DATA_FILE, REF_FILE_SAMPLE_DATA);
        RefSnvCountsFile = getRefDataFile(configBuilder, REF_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS);
        RefSigContributionFile = getRefDataFile(configBuilder, REF_SIG_CONTRIB_FILE, REF_FILE_SIG_PERC);
        RefFeaturePrevFile = getRefDataFile(configBuilder, REF_FEAT_PREV_FILE, REF_FILE_FEATURE_PREV);
        RefTraitPercFile = getRefDataFile(configBuilder, REF_TRAIT_PERC_FILE, REF_FILE_TRAIT_PERC);
        RefTraitRateFile = getRefDataFile(configBuilder, REF_TRAIT_RATE_FILE, REF_FILE_TRAIT_RATES);
        RefGenderRateFile = getRefDataFile(configBuilder, REF_GENDER_RATE_FILE, REF_FILE_GENDER_RATES);
        RefSvPercFile = getRefDataFile(configBuilder, REF_SV_PERC_FILE, REF_FILE_SV_PERC);
        RefSnvCancerPosFreqFile = getRefDataFile(configBuilder, REF_SNV_CANCER_POS_FREQ_FILE, REF_FILE_CANCER_POS_FREQ_COUNTS);
        RefSnvSamplePosFreqFile = getRefDataFile(configBuilder, REF_SNV_SAMPLE_POS_FREQ_FILE, REF_FILE_SAMPLE_POS_FREQ_COUNTS, true);
        RefDriverAvgFile = getRefDataFile(configBuilder, REF_DRIVER_AVG_FILE, REF_FILE_DRIVER_AVG);
        RefSnvSignaturesFile = getRefDataFile(configBuilder, REF_SNV_SIGNATURES_FILE, REF_FILE_SNV_SIGNATURES);

        RefGeneExpCancerFile = getRefDataFile(configBuilder, REF_RNA_GENE_EXP_CANCER_FILE, REF_FILE_GENE_EXP_CANCER, true);
        RefGeneExpSampleFile = getRefDataFile(configBuilder, REF_RNA_GENE_EXP_SAMPLE_FILE, REF_FILE_GENE_EXP_SAMPLE, true);
        RefAltSjCancerFile = getRefDataFile(configBuilder, REF_RNA_ALT_SJ_CANCER_FILE, REF_FILE_ALT_SJ_CANCER, true);
        RefAltSjSampleFile = getRefDataFile(configBuilder, REF_RNA_ALT_SJ_SAMPLE_FILE, REF_FILE_ALT_SJ_SAMPLE, true);

        TestRefData = configBuilder.hasFlag(TEST_REF_SAMPLE_DATA);

        // use cases for loading sample data:
        // 1. DB - sourced
        // 2. single sample - uses pipeline names for each type (eg Linx, Purple, Isofox)
        // 3. Cohort and Reference data files - if 'use_ref_sample_data' specified, then run Cuppa over ref & cohort files
        // 4. Non-reference cohort files specified manually
        // 5. Sample data not supplied for a given category

        if(TestRefData)
        {
            SampleDataFile = RefSampleDataFile;

            RefSampleTraitsFile = getRefDataFile(configBuilder, REF_COHORT_SAMPLE_TRAITS_FILE, COHORT_REF_TRAITS_DATA_FILE);
            RefSampleFeatureFile = getRefDataFile(configBuilder, REF_COHORT_FEATURES_FILE, COHORT_REF_FEATURE_DATA_FILE);
            RefSampleSigContribFile = getRefDataFile(configBuilder, REF_COHORT_SIG_CONTRIBS_FILE, COHORT_REF_SIG_DATA_FILE);
            RefSampleSvFile = getRefDataFile(configBuilder, REF_COHORT_SV_DATA_FILE, COHORT_REF_SV_DATA_FILE);

            if(configBuilder.hasValue(SAMPLE))
            {
                CUP_LOGGER.info("testing single reference sample({})", configBuilder.getValue(SAMPLE));
            }
            else
            {
                CUP_LOGGER.info("testing all reference samples");
            }

            SampleDataDir = "";
            LinxDir = "";
            PurpleDir = "";
            IsofoxDir = "";
            VirusDir = "";
            DbAccess = null;
        }
        else
        {
            SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG, ""));
            SampleDataFile = configBuilder.getValue(SAMPLE_DATA_FILE, "");

            LinxDir = checkAddDirSeparator(configBuilder.getValue(LINX_DIR_CFG, ""));
            PurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG, ""));
            IsofoxDir = checkAddDirSeparator(configBuilder.getValue(ISOFOX_DIR_CFG, ""));
            VirusDir = checkAddDirSeparator(configBuilder.getValue(VIRUS_DIR_CFG, ""));

            if(configBuilder.hasValue(SAMPLE))
            {
                CUP_LOGGER.info("testing single sample({})", configBuilder.getValue(SAMPLE));
            }
            else if(configBuilder.hasValue(SAMPLE_DATA_DIR_CFG))
            {
                CUP_LOGGER.info("testing samples from file: {}", SampleDataFile);
            }
            else
            {
                CUP_LOGGER.error("missing {}, non-ref cohort {} or {} config", SAMPLE, SAMPLE_DATA_FILE, TEST_REF_SAMPLE_DATA);
            }

            RefSampleTraitsFile = "";
            RefSampleSigContribFile = "";
            RefSampleFeatureFile = "";
            RefSampleSvFile = "";

            DbAccess = createDatabaseAccess(configBuilder);
        }

        NoiseAdjustments = new NoiseRefCache(RefDataDir);
        NoiseAdjustments.loadNoiseAllocations(configBuilder.getValue(NOISE_ALLOCATIONS));

        FeatureDampenFactor = configBuilder.getDecimal(FEATURE_DAMPEN_FACTOR);
        CombinedDampenFactor = configBuilder.getDecimal(COMBINED_DAMPEN_FACTOR);

        NoSubtypeCollapse = configBuilder.hasFlag(NO_SUBTYPE_COLLAPSE);

        OutputDir = parseOutputDir(configBuilder);
        OutputFileId = configBuilder.getValue(OUTPUT_ID, "");
        Threads = parseThreads(configBuilder);

        WriteSimilarities = configBuilder.hasFlag(WRITE_SIMS);
        WriteCondensed = configBuilder.hasFlag(WRITE_CONDENSED);
        WriteDetailedScores = configBuilder.hasFlag(WRITE_DETAILED_SCORES);
        CreatePdf = configBuilder.hasFlag(CREATE_PDF);
    }

    private String getRefDataFile(final ConfigBuilder configBuilder, final String configStr, final String defaultFilename)
    {
        return getRefDataFile(configBuilder, configStr, defaultFilename, false);
    }

    private String getRefDataFile(final ConfigBuilder configBuilder, final String configStr, final String defaultFilename, boolean checkZipped)
    {
        String refFilename = configBuilder.hasValue(configStr) ? configBuilder.getValue(configStr) : RefDataDir + defaultFilename;

        if(checkZipped && !Files.exists(Paths.get(refFilename)) && Files.exists(Paths.get(refFilename + ".gz")))
            return refFilename + ".gz";

        return refFilename;
    }

    public boolean isValid()
    {
        return !OutputDir.isEmpty();
    }

    public boolean runClassifier(final CategoryType type) { return classifierEnabled(type, Categories); }

    public static boolean classifierEnabled(final CategoryType type, final List<CategoryType> categories)
    {
        return categories.isEmpty() || categories.contains(type);
    }

    public String formOutputFilename(final String fileId)
    {
        String outputFile = OutputDir + "CUP";

        if(!OutputFileId.isEmpty())
            outputFile += "." + OutputFileId;

        return outputFile + "." + fileId + ".csv";
    }

    public String getLinxDataDir(final String sampleId) { return getSampleDataDir(sampleId, LinxDir); }
    public String getPurpleDataDir(final String sampleId) { return getSampleDataDir(sampleId, PurpleDir); }
    public String getVirusDataDir(final String sampleId) { return getSampleDataDir(sampleId, VirusDir); }
    public String getIsofoxDataDir(final String sampleId) { return getSampleDataDir(sampleId, IsofoxDir); }

    private String  getSampleDataDir(final String sampleId, final String specificDir)
    {
        return !SampleDataDir.isEmpty() ? formSamplePath(SampleDataDir, sampleId) : formSamplePath(specificDir, sampleId);
    }

    public static String formSamplePath(final String samplePath, final String sampleId)
    {
        return convertWildcardSamplePath(samplePath, sampleId);
    }

    public static List<CategoryType> configCategories(final ConfigBuilder configBuilder)
    {
        List<CategoryType> categories = Lists.newArrayList();

        if(configBuilder.hasValue(CATEGORIES))
        {
            if(configBuilder.getValue(CATEGORIES).equals(ALL_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> x != COMBINED).forEach(x -> categories.add(x));
            }
            else if(configBuilder.getValue(CATEGORIES).equals(DNA_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> isDna(x)).forEach(x -> categories.add(x));
            }
            else if(configBuilder.getValue(CATEGORIES).equals(RNA_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> isRna(x)).forEach(x -> categories.add(x));
            }
            else
            {
                final String[] categoryStrings = configBuilder.getValue(CATEGORIES).split(SUBSET_DELIM);
                Arrays.stream(categoryStrings).forEach(x -> categories.add(CategoryType.valueOf(x)));
            }
        }
        else
        {
            // just DNA by default
            Arrays.stream(CategoryType.values()).filter(x -> isDna(x)).forEach(x -> categories.add(x));
        }

        return categories;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        StringJoiner categories = new StringJoiner(",");
        Arrays.stream(CategoryType.values()).filter(x -> x != COMBINED).forEach(x -> categories.add(x.toString()));

        configBuilder.addConfigItem(
                CATEGORIES, false,
                format("Categories for analysis: %s, %s, %s or sub-group from [%s]",
                        ALL_CATEGORIES, DNA_CATEGORIES, RNA_CATEGORIES, categories));

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addInteger(SAMPLE_RNA_LENGTH, "Specific sample RNA read length", DEFAULT_RNA_LENGTH);
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);
        addPipelineDirectories(configBuilder);

        configBuilder.addPath(SAMPLE_DATA_FILE, false, "Sample data file");

        configBuilder.addPath(REF_COHORT_SIG_CONTRIBS_FILE, false, "Cohort ref sample signature contributions");
        configBuilder.addPath(REF_COHORT_FEATURES_FILE, false, "Cohort ref features file (drivers, fusions and viruses)");
        configBuilder.addPath(REF_COHORT_SAMPLE_TRAITS_FILE, false, "Cohort sample traits file");
        configBuilder.addPath(REF_COHORT_SV_DATA_FILE, false, "Cohort SV data");

        configBuilder.addPath(REF_DATA_DIR, true, "Reference data directory");
        configBuilder.addFlag(TEST_REF_SAMPLE_DATA, "In cohort-mode, run Cuppa using all ref sample data files");

        configBuilder.addPath(REF_SAMPLE_DATA_FILE, false, "Reference sample data, default: " + REF_FILE_SAMPLE_DATA);
        configBuilder.addPath(REF_SNV_COUNTS_FILE, false, "Reference SNV sample counts, default: " + REF_FILE_SNV_COUNTS);
        configBuilder.addPath(REF_SIG_CONTRIB_FILE, false, "SNV signatures, default: " + REF_FILE_SIG_PERC);
        configBuilder.addPath(REF_FEAT_PREV_FILE, false, "Reference driver prevalence, default: " + REF_FILE_FEATURE_PREV);
        configBuilder.addPath(REF_SV_PERC_FILE, false, "Reference SV percentiles file, default: " + REF_FILE_SV_PERC);
        configBuilder.addPath(REF_TRAIT_PERC_FILE, false, "Reference traits percentiles file, default: " + REF_FILE_TRAIT_PERC);
        configBuilder.addPath(REF_TRAIT_RATE_FILE, false, "Reference traits rates file, default: " + REF_FILE_TRAIT_RATES);
        configBuilder.addPath(REF_GENDER_RATE_FILE, false, "Reference gender rates file, default: " + REF_FILE_GENDER_RATES);
        configBuilder.addPath(REF_SNV_CANCER_POS_FREQ_FILE, false, "Reference SNV cancer position frequency file, default: " + REF_FILE_CANCER_POS_FREQ_COUNTS);
        configBuilder.addPath(REF_SNV_SAMPLE_POS_FREQ_FILE, false, "Reference SNV sample position frequency file, default: " + REF_FILE_SAMPLE_POS_FREQ_COUNTS);
        configBuilder.addPath(REF_DRIVER_AVG_FILE, false, "Reference features per sample file, default: " + REF_FILE_DRIVER_AVG);
        configBuilder.addPath(REF_SNV_SIGNATURES_FILE, false, "Reference SNV signatures, default: " + REF_FILE_SNV_SIGNATURES);
        configBuilder.addPath(REF_RNA_GENE_EXP_CANCER_FILE, false, "Reference RNA cancer gene expression file, default: " + REF_FILE_GENE_EXP_CANCER);
        configBuilder.addPath(REF_RNA_GENE_EXP_SAMPLE_FILE, false, "Reference RNA sample gene expression file, default: " + REF_FILE_GENE_EXP_SAMPLE);
        configBuilder.addPath(REF_RNA_ALT_SJ_CANCER_FILE, false, "Reference RNA alternative splice-junction cancer file, default: " + REF_FILE_ALT_SJ_CANCER);
        configBuilder.addPath(REF_RNA_ALT_SJ_SAMPLE_FILE, false, "Reference RNA alternative splice-junction sample file, default: " + REF_FILE_ALT_SJ_SAMPLE);
        configBuilder.addConfigItem(NOISE_ALLOCATIONS, NOISE_ALLOCATIONS_DESC);
        configBuilder.addFlag(NO_SUBTYPE_COLLAPSE, "Keep cancer sub-types separated in final classifiers");

        configBuilder.addDecimal(COMBINED_DAMPEN_FACTOR,"Combined classifier dampening factor", COMBINED_DAMPEN_FACTOR_DEFAULT);

        configBuilder.addFlag(WRITE_SIMS, "Write top-20 CSS similarities to file");
        configBuilder.addFlag(WRITE_DETAILED_SCORES, "Cohort-only - write detailed (non-classifier) data");
        configBuilder.addFlag(WRITE_CONDENSED, "Write sample results as single line");
        configBuilder.addFlag(CREATE_PDF, "Call CUP Report ccript to generate PDF");

        addDatabaseCmdLineArgs(configBuilder, false);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        GeneExpressionClassifier.addCmdLineArgs(configBuilder);
        AltSjClassifier.addCmdLineArgs(configBuilder);
        SomaticClassifier.registerConfig(configBuilder);
        FeatureClassifier.registerConfig(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }

    public CuppaConfig()
    {
        RefGenVersion = V37;
        Categories = Lists.newArrayList();
        RefDataDir = "";
        RefSampleDataFile = "";
        RefSnvCountsFile = "";
        RefSvPercFile = "";
        RefSigContributionFile = "";
        RefFeaturePrevFile = "";
        RefTraitPercFile = "";
        RefTraitRateFile = "";
        RefGenderRateFile = "";
        RefSnvCancerPosFreqFile = "";
        RefSnvSamplePosFreqFile = "";
        RefDriverAvgFile = "";
        RefSnvSignaturesFile = "";

        RefGeneExpCancerFile = "";
        RefGeneExpSampleFile = "";
        RefAltSjCancerFile = "";
        RefAltSjSampleFile = "";

        TestRefData = false;

        // sample data, if not sourced from the database
        SampleDataDir = "";
        LinxDir = "";
        PurpleDir = "";
        IsofoxDir = "";
        VirusDir = "";

        SampleDataFile = "";
        RefSampleFeatureFile = "";
        RefSampleTraitsFile = "";
        RefSampleSigContribFile = "";
        RefSampleSvFile = "";

        NoiseAdjustments = new NoiseRefCache(null);
        NoSubtypeCollapse = false;

        FeatureDampenFactor = FEATURE_DAMPEN_FACTOR_DEFAULT;
        CombinedDampenFactor = COMBINED_DAMPEN_FACTOR_DEFAULT;

        DbAccess = null;
        WriteSimilarities = false;
        WriteDetailedScores = false;
        WriteCondensed = false;
        CreatePdf = false;
        OutputDir = "";
        OutputFileId = "";
        Threads = 0;
    }
}
