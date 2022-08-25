package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.ALL_CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.NOISE_ALLOCATIONS;
import static com.hartwig.hmftools.cup.CuppaConfig.NOISE_ALLOCATIONS_DESC;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_COUNTS_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_SAMPLE_POS_FREQ_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.configCategories;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.feature.RefFeatures;
import com.hartwig.hmftools.cup.rna.RefGeneExpression;
import com.hartwig.hmftools.cup.somatics.RefSomatics;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RefDataConfig
{
    public final List<CategoryType> Categories;

    public final String OutputDir;

    // reference data
    public final String SampleDataFile;
    public final String CohortSampleTraitsFile;
    public final String CohortSampleSvDataFile;
    public final String CohortSigContribsFile;
    public final String CohortFeaturesFile;
    public final String GenPosMatrixFile;
    public final String Snv96MatrixFile;
    public final String GeneExpMatrixFile;
    public final String AltSjMatrixFile;
    public final String FeatureOverrideFile;

    // pipeline directories, accepting wildcards
    public final String LinxDir;
    public final String PurpleDir;
    public final String IsofoxDir;

    public final DatabaseAccess DbAccess;

    public final NoiseRefCache NoiseAdjustments;

    public final boolean WriteCohortFiles; // re-write data sourced from database or flat files into single cohort files

    // config strings

    // cohort files instead of extracting data from database
    public static final String REF_COHORT_SAMPLE_TRAITS_FILE = "cohort_sample_traits_file";
    public static final String REF_COHORT_SIG_CONTRIBS_FILE = "cohort_sig_contribs_file";
    public static final String REF_COHORT_SV_DATA_FILE = "cohort_sv_data_file";
    public static final String REF_COHORT_FEATURES_FILE = "cohort_features_file";

    private static final String REF_GENE_EXP_DATA_FILE = "ref_gene_exp_file";
    private static final String REF_ALT_SJ_DATA_FILE = "ref_alt_sj_file";

    // pipeline files from Linx, Isofox and Purple
    public static final String LINX_DIR = "linx_dir";
    public static final String PURPLE_DIR = "purple_dir";
    public static final String ISOFOX_DIR = "isofox_dir";

    private static final String REF_FEATURE_OVERRIDE_FILE = "feature_override_file";
    public static final String GENDER_RATES = "gender_rates";
    private static final String WRITE_COHORT_FILES = "write_cohort_files";

    private static final String FILE_DELIM = ",";

    public RefDataConfig(final CommandLine cmd)
    {
        Categories = configCategories(cmd);

        CUP_LOGGER.info("build ref data for classifiers: {}", Categories.isEmpty() ? ALL_CATEGORIES : Categories.toString());

        SampleDataFile = cmd.getOptionValue(REF_SAMPLE_DATA_FILE, "");
        CohortSampleTraitsFile = cmd.getOptionValue(REF_COHORT_SAMPLE_TRAITS_FILE, "");
        CohortSigContribsFile = cmd.getOptionValue(REF_COHORT_SIG_CONTRIBS_FILE, "");
        CohortSampleSvDataFile = cmd.getOptionValue(REF_COHORT_SV_DATA_FILE, "");
        CohortFeaturesFile = cmd.getOptionValue(REF_COHORT_FEATURES_FILE, "");

        LinxDir = cmd.getOptionValue(LINX_DIR, "");
        PurpleDir = cmd.getOptionValue(PURPLE_DIR, "");
        IsofoxDir = cmd.getOptionValue(ISOFOX_DIR, "");

        GenPosMatrixFile = cmd.getOptionValue(REF_SNV_SAMPLE_POS_FREQ_FILE, "");
        Snv96MatrixFile = cmd.getOptionValue(REF_SNV_COUNTS_FILE, "");

        GeneExpMatrixFile = cmd.getOptionValue(REF_GENE_EXP_DATA_FILE, "");
        AltSjMatrixFile = cmd.getOptionValue(REF_ALT_SJ_DATA_FILE, "");

        FeatureOverrideFile = cmd.getOptionValue(REF_FEATURE_OVERRIDE_FILE, "");

        DbAccess = createDatabaseAccess(cmd);

        OutputDir = parseOutputDir(cmd);

        NoiseAdjustments = new NoiseRefCache(OutputDir);
        NoiseAdjustments.loadNoiseAllocations(cmd.getOptionValue(NOISE_ALLOCATIONS));

        WriteCohortFiles = cmd.hasOption(WRITE_COHORT_FILES);
    }

    public static final List<String> parseFileSet(final String filenames)
    {
        return Arrays.stream(filenames.split(FILE_DELIM, -1)).collect(Collectors.toList());
    }

    public static void addPipelineDirectories(final Options options)
    {
        options.addOption(LINX_DIR, true, "Path to Linx data files for sample, wildcards allowed");
        options.addOption(PURPLE_DIR, true, "Path to Purple data files for sample, wildcards allowed");
        options.addOption(ISOFOX_DIR, true, "Path to Purple data files for sample, wildcards allowed");
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CATEGORIES, true, "Categories to build ref data for");

        options.addOption(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");
        options.addOption(REF_SNV_SAMPLE_POS_FREQ_FILE, true, "Ref SNV position frequency matrix data file");
        options.addOption(REF_SNV_COUNTS_FILE, true, "Ref SNV trinucleotide matrix data file");

        options.addOption(REF_COHORT_FEATURES_FILE, true, "Ref sample features data file");
        options.addOption(REF_COHORT_SAMPLE_TRAITS_FILE, true, "Ref sample cohort traits file");
        options.addOption(REF_COHORT_SIG_CONTRIBS_FILE, true, "Ref sample cohort signature contributions file");
        options.addOption(REF_COHORT_SV_DATA_FILE, true, "Ref sample cohort SV file");

        addPipelineDirectories(options);

        options.addOption(REF_GENE_EXP_DATA_FILE, true, "Ref sample RNA gene expression cohort data file");
        options.addOption(REF_ALT_SJ_DATA_FILE, true, "Ref sample RNA alternate splice junction cohort data file");

        options.addOption(REF_FEATURE_OVERRIDE_FILE, true, "Ref feature override data file");

        options.addOption(NOISE_ALLOCATIONS, true, NOISE_ALLOCATIONS_DESC);
        options.addOption(GENDER_RATES, true, "Gender-rate overrides - format CancerType;MalePerc;FemalePerc, etc");
        options.addOption(WRITE_COHORT_FILES, false, "Re-write ref data as cohort files");

        addDatabaseCmdLineArgs(options);

        addLoggingOptions(options);
        addOutputOptions(options);

        RefGeneExpression.addCmdLineArgs(options);
        RefSomatics.addCmdLineArgs(options);
        RefFeatures.addCmdLineArgs(options);
    }

}
