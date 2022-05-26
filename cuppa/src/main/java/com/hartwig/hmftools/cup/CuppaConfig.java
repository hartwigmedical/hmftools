package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_FEATURE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_SIG_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_SV_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.COHORT_REF_FILE_TRAITS_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_CANCER;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_CANCER_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_COPY_NUMBER_PROFILE;
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
import static com.hartwig.hmftools.cup.common.CategoryType.ALL_CATEGORIES;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.COMBINED;
import static com.hartwig.hmftools.cup.common.CategoryType.DNA_CATEGORIES;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.RNA_CATEGORIES;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.CategoryType.isDna;
import static com.hartwig.hmftools.cup.common.CategoryType.isRna;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.ISOFOX_DIR;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.LINX_DIR;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.PURPLE_DIR;
import static com.hartwig.hmftools.cup.ref.RefDataConfig.addPipelineDirectories;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.rna.AltSjClassifier;
import com.hartwig.hmftools.cup.rna.GeneExpressionClassifier;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;
import com.hartwig.hmftools.cup.traits.SampleTraitClassifier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CuppaConfig
{
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
    public final String IsofoxDir;

    // cohort files, formed during ref data building
    public final String SampleDataFile;
    public final String SampleFeatureFile;
    public final String SampleTraitsFile;
    public final String SampleSnvCountsFile;
    public final String SampleSnvPosFreqFile;
    public final String SampleSigContribFile;
    public final String SampleSvFile;
    public final String SampleGeneExpFile;
    public final String SampleAltSjFile;

    public final NoiseRefCache NoiseAdjustments;

    // database access
    public final DatabaseAccess DbAccess;

    public final boolean WriteSimilarities;
    public final boolean WriteDetailedScores;
    public final boolean WriteCondensed;

    public final String OutputDir;
    public final String OutputFileId;
    public final int Threads;

    // config strings
    public static final String CATEGORIES = "categories";

    public static final String SAMPLE_DATA_DIR = "sample_data_dir";

    public static final String SPECIFIC_SAMPLE_DATA = "sample_data";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    private static final String SAMPLE_FEAT_FILE = "sample_feature_file";
    private static final String SAMPLE_TRAITS_FILE = "sample_traits_file";
    private static final String SAMPLE_SNV_COUNTS_FILE = "sample_snv_counts_file";
    private static final String SAMPLE_SNV_POS_FREQ_FILE = "sample_snv_pos_freq_file";
    private static final String SAMPLE_SIG_CONTRIB_FILE = "sample_sig_contrib_file";
    private static final String SAMPLE_SOMATIC_VCF = "sample_somatic_vcf";
    private static final String SAMPLE_SV_FILE = "sample_sv_file";
    private static final String SAMPLE_GENE_EXP_FILE = "sample_gene_exp_file";
    private static final String SAMPLE_ALT_SJ_FILE = "sample_alt_sj_matrix_file";

    public static final String REF_DATA_DIR = "ref_data_dir";
    public static final String USE_REF_SAMPLE_DATA = "use_ref_sample_data";

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

    public static final String WRITE_SIMS = "write_similarities";
    public static final String WRITE_DETAILED_SCORES = "write_detailed_scores";
    public static final String WRITE_CONDENSED = "write_condensed";

    public static final String OUTPUT_FILE_ID = "output_id";
    public static final String LOG_DEBUG = "log_debug";
    public static final String THREADS = "threads";

    public static final Logger CUP_LOGGER = LogManager.getLogger(CuppaConfig.class);

    // file fields
    public static final String FLD_SAMPLE_ID = "SampleId";
    public static final String FLD_CANCER_TYPE = "CancerType";
    public static final String FLD_CANCER_SUBTYPE = "CancerSubtype";
    public static final String FLD_RNA_READ_LENGTH = "RnaReadLength";
    public static final String CANCER_SUBTYPE_OTHER = "Other";
    public static final String DATA_DELIM = ",";
    public static final String SUBSET_DELIM = ";";

    public CuppaConfig(final CommandLine cmd)
    {
        Categories = configCategories(cmd);

        CUP_LOGGER.info("running classifiers: {}", Categories.isEmpty() ? ALL_CATEGORIES : Categories.toString());

        RefDataDir = checkAddDirSeparator(cmd.getOptionValue(REF_DATA_DIR, ""));

        RefSampleDataFile = getRefDataFile(cmd, REF_SAMPLE_DATA_FILE, REF_FILE_SAMPLE_DATA);
        RefSnvCountsFile = getRefDataFile(cmd, REF_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS);
        RefSigContributionFile = getRefDataFile(cmd, REF_SIG_CONTRIB_FILE, REF_FILE_SIG_PERC);
        RefFeaturePrevFile = getRefDataFile(cmd, REF_FEAT_PREV_FILE, REF_FILE_FEATURE_PREV);
        RefTraitPercFile = getRefDataFile(cmd, REF_TRAIT_PERC_FILE, REF_FILE_TRAIT_PERC);
        RefTraitRateFile = getRefDataFile(cmd, REF_TRAIT_RATE_FILE, REF_FILE_TRAIT_RATES);
        RefGenderRateFile = getRefDataFile(cmd, REF_GENDER_RATE_FILE, REF_FILE_GENDER_RATES);
        RefSvPercFile = getRefDataFile(cmd, REF_SV_PERC_FILE, REF_FILE_SV_PERC);
        RefSnvCancerPosFreqFile = getRefDataFile(cmd, REF_SNV_CANCER_POS_FREQ_FILE, REF_FILE_CANCER_POS_FREQ_COUNTS);
        RefSnvSamplePosFreqFile = getRefDataFile(cmd, REF_SNV_SAMPLE_POS_FREQ_FILE, REF_FILE_SAMPLE_POS_FREQ_COUNTS);
        RefDriverAvgFile = getRefDataFile(cmd, REF_DRIVER_AVG_FILE, REF_FILE_DRIVER_AVG);
        RefSnvSignaturesFile = getRefDataFile(cmd, REF_SNV_SIGNATURES_FILE, REF_FILE_SNV_SIGNATURES);

        RefGeneExpCancerFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_CANCER_FILE, REF_FILE_GENE_EXP_CANCER);
        RefGeneExpSampleFile = getRefDataFile(cmd, REF_RNA_GENE_EXP_SAMPLE_FILE, REF_FILE_GENE_EXP_SAMPLE);
        RefAltSjCancerFile = getRefDataFile(cmd, REF_RNA_ALT_SJ_CANCER_FILE, REF_FILE_ALT_SJ_CANCER);
        RefAltSjSampleFile = getRefDataFile(cmd, REF_RNA_ALT_SJ_SAMPLE_FILE, REF_FILE_ALT_SJ_SAMPLE);

        boolean useRefData = cmd.hasOption(USE_REF_SAMPLE_DATA);
        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR, ""));
        LinxDir = checkAddDirSeparator(cmd.getOptionValue(LINX_DIR, ""));
        PurpleDir = checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR, ""));
        IsofoxDir = checkAddDirSeparator(cmd.getOptionValue(ISOFOX_DIR, ""));

        // use cases for loading sample data:
        // 1. DB - sourced
        // 2. single sample - uses pipeline names for each type (eg Linx, Purple, Isofox)
        // 3. Cohort and Reference data files - if 'use_ref_sample_data' specified, then run Cuppa over ref & cohort files
        // 4. Cohort files specified manually
        // 5. Sample data not supplied for a given data type

        SampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE, "");

        // for cohort mode, then re-testing ref data, the cohort files are used here:
        SampleTraitsFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_TRAITS_FILE, COHORT_REF_FILE_TRAITS_DATA_FILE, SAMPLE_TRAIT);
        SampleFeatureFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_FEAT_FILE, COHORT_REF_FILE_FEATURE_DATA_FILE, FEATURE);
        SampleSigContribFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_SIG_CONTRIB_FILE, COHORT_REF_FILE_SIG_DATA_FILE, SNV);
        SampleSvFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_SV_FILE, COHORT_REF_FILE_SV_DATA_FILE, SV);

        SampleSnvCountsFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_SNV_COUNTS_FILE, REF_FILE_SNV_COUNTS, SNV);
        SampleSnvPosFreqFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_SNV_POS_FREQ_FILE, REF_FILE_SAMPLE_POS_FREQ_COUNTS, SNV);
        SampleGeneExpFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_GENE_EXP_FILE, REF_FILE_GENE_EXP_SAMPLE, GENE_EXP);
        SampleAltSjFile = getCohortSampleDataFile(cmd, useRefData, SAMPLE_ALT_SJ_FILE, REF_FILE_ALT_SJ_SAMPLE, ALT_SJ);

        NoiseAdjustments = new NoiseRefCache(RefDataDir);
        NoiseAdjustments.loadNoiseAllocations(cmd.getOptionValue(NOISE_ALLOCATIONS));

        OutputDir = parseOutputDir(cmd);
        OutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID, "");
        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));

        WriteSimilarities = cmd.hasOption(WRITE_SIMS);
        WriteCondensed = cmd.hasOption(WRITE_CONDENSED);
        WriteDetailedScores = cmd.hasOption(WRITE_DETAILED_SCORES);

        DbAccess = createDatabaseAccess(cmd);
    }

    private String getRefDataFile(final CommandLine cmd, final String configStr, final String defaultFilename)
    {
        if(cmd.hasOption(configStr))
            return cmd.getOptionValue(configStr);

        return RefDataDir + defaultFilename;
    }

    private String getCohortSampleDataFile(
            final CommandLine cmd, boolean useRefDataFile, final String configStr, final String defaultFilename, final CategoryType category)
    {
        if(cmd.hasOption(SAMPLE_DATA_DIR) && cmd.hasOption(SPECIFIC_SAMPLE_DATA))
            return "";

        if(cmd.hasOption(configStr))
            return cmd.getOptionValue(configStr);

        if(!useRefDataFile)
            return ""; // meaning this data type is not loaded

        if(!runClassifier(category)) // ignores ref data file name since classifier won't be run
            return "";

        return defaultFilename.startsWith(RefDataDir) ? defaultFilename : RefDataDir + defaultFilename;
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
    public String getIsofoxDataDir(final String sampleId) { return getSampleDataDir(sampleId, IsofoxDir); }

    private String getSampleDataDir(final String sampleId, final String specificDir)
    {
        return !SampleDataDir.isEmpty() ? formSamplePath(SampleDataDir, sampleId) : formSamplePath(specificDir, sampleId);
    }

    public static String formSamplePath(final String samplePath, final String sampleId)
    {
        if(!samplePath.contains("*"))
            return samplePath;

        return samplePath.replaceAll("\\*", sampleId);
    }

    public static List<CategoryType> configCategories(final CommandLine cmd)
    {
        List<CategoryType> categories = Lists.newArrayList();

        if(cmd.hasOption(CATEGORIES))
        {
            if(cmd.getOptionValue(CATEGORIES).equals(ALL_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> x != COMBINED).forEach(x -> categories.add(x));
            }
            else if(cmd.getOptionValue(CATEGORIES).equals(DNA_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> isDna(x)).forEach(x -> categories.add(x));
            }
            else if(cmd.getOptionValue(CATEGORIES).equals(RNA_CATEGORIES))
            {
                Arrays.stream(CategoryType.values()).filter(x -> isRna(x)).forEach(x -> categories.add(x));
            }
            else
            {
                final String[] categoryStrings = cmd.getOptionValue(CATEGORIES).split(SUBSET_DELIM);
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

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(
                CATEGORIES, true,
                "Categories to run analysis on, separated by ';' from: SNV, SV, SAMPLE_TRAIT, GENE_EXP and FEATURE");

        options.addOption(SPECIFIC_SAMPLE_DATA, true, "Specific sample in form 'SampleId;CancerType;CancerSubtype' (last 2 optional)");
        options.addOption(SAMPLE_DATA_DIR, true, "Directory containing standard sample files from pipeline");
        addPipelineDirectories(options);

        options.addOption(SAMPLE_DATA_FILE, true, "Sample data file");

        options.addOption(SAMPLE_SNV_COUNTS_FILE, true, "Sample SNV counts");
        options.addOption(SAMPLE_SNV_POS_FREQ_FILE, true, "Sample SNV position frequence counts");
        options.addOption(SAMPLE_SIG_CONTRIB_FILE, true, "Sample signature contributions");

        options.addOption(SAMPLE_FEAT_FILE, true, "Cohort features file (drivers, fusions and viruses)");
        options.addOption(SAMPLE_TRAITS_FILE, true, "Cohort sample traits file");
        options.addOption(SAMPLE_SV_FILE, true, "Cohort SV data");
        options.addOption(SAMPLE_GENE_EXP_FILE, true, "Cohort sample RNA gene expression TPMs");
        options.addOption(SAMPLE_ALT_SJ_FILE, true, "Cohort sample RNA alt-SJ frag counts matrix file");
        options.addOption(SAMPLE_SOMATIC_VCF, true, "Sample somatic VCF");

        options.addOption(REF_DATA_DIR, true, "Reference data directory");
        options.addOption(USE_REF_SAMPLE_DATA, false, "In cohort-mode, run Cuppa using all ref sample data files");

        options.addOption(REF_SAMPLE_DATA_FILE, true, "Reference sample data, default: " + REF_FILE_SAMPLE_DATA);
        options.addOption(REF_SNV_COUNTS_FILE, true, "Reference SNV sample counts, default: " + REF_FILE_SNV_COUNTS);
        options.addOption(REF_SIG_CONTRIB_FILE, true, "SNV signatures, default: " + REF_FILE_SIG_PERC);
        options.addOption(REF_FEAT_PREV_FILE, true, "Reference driver prevalence, default: " + REF_FILE_FEATURE_PREV);
        options.addOption(REF_SV_PERC_FILE, true, "Reference SV percentiles file, default: " + REF_FILE_SV_PERC);
        options.addOption(REF_TRAIT_PERC_FILE, true, "Reference traits percentiles file, default: " + REF_FILE_TRAIT_PERC);
        options.addOption(REF_TRAIT_RATE_FILE, true, "Reference traits rates file, default: " + REF_FILE_TRAIT_RATES);
        options.addOption(REF_SNV_CANCER_POS_FREQ_FILE, true, "Reference SNV cancer position frequency file, default: " + REF_FILE_CANCER_POS_FREQ_COUNTS);
        options.addOption(REF_SNV_SAMPLE_POS_FREQ_FILE, true, "Reference SNV sample position frequency file, default: " + REF_FILE_SAMPLE_POS_FREQ_COUNTS);
        options.addOption(REF_DRIVER_AVG_FILE, true, "Reference features per sample file, default: " + REF_FILE_DRIVER_AVG);
        options.addOption(REF_SNV_SIGNATURES_FILE, true, "Reference SNV signatures, default: " + REF_FILE_SNV_SIGNATURES);
        options.addOption(REF_RNA_GENE_EXP_CANCER_FILE, true, "Reference RNA cancer gene expression file, default: " + REF_FILE_GENE_EXP_CANCER);
        options.addOption(REF_RNA_GENE_EXP_SAMPLE_FILE, true, "Reference RNA sample gene expression file, default: " + REF_FILE_GENE_EXP_SAMPLE);
        options.addOption(REF_RNA_ALT_SJ_CANCER_FILE, true, "Reference RNA alternative splice-junction cancer file, default: " + REF_FILE_ALT_SJ_CANCER);
        options.addOption(REF_RNA_ALT_SJ_SAMPLE_FILE, true, "Reference RNA alternative splice-junction sample file, default: " + REF_FILE_ALT_SJ_SAMPLE);
        options.addOption(NOISE_ALLOCATIONS, true, NOISE_ALLOCATIONS_DESC);

        options.addOption(WRITE_SIMS, false, "Write top-20 CSS similarities to file");
        options.addOption(WRITE_DETAILED_SCORES, false, "Cohort-only - write detailed (non-classifier) data");
        options.addOption(WRITE_CONDENSED, false, "Write sample results as single line");

        addDatabaseCmdLineArgs(options);
        GeneExpressionClassifier.addCmdLineArgs(options);
        AltSjClassifier.addCmdLineArgs(options);
        SomaticClassifier.addCmdLineArgs(options);
        FeatureClassifier.addCmdLineArgs(options);
        SampleTraitClassifier.addCmdLineArgs(options);

        options.addOption(OUTPUT_DIR, true, "Path to output files");
        options.addOption(OUTPUT_FILE_ID, true, "Output file ID");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(THREADS, true, "Number of threads");
    }

    public CuppaConfig()
    {
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

        // sample data, if not sourced from the database
        SampleDataDir = "";
        LinxDir = "";
        PurpleDir = "";
        IsofoxDir = "";

        SampleDataFile = "";
        SampleFeatureFile = "";
        SampleTraitsFile = "";
        SampleSnvCountsFile = "";
        SampleSnvPosFreqFile = "";
        SampleSigContribFile = "";
        SampleSvFile = "";
        SampleGeneExpFile = "";
        SampleAltSjFile = "";

        NoiseAdjustments = new NoiseRefCache(null);

        DbAccess = null;
        WriteSimilarities = false;
        WriteDetailedScores = false;
        WriteCondensed = false;
        OutputDir = "";
        OutputFileId = "";
        Threads = 0;
    }
}
