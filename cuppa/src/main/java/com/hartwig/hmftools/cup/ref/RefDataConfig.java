package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.ALL_CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.NOISE_ALLOCATIONS;
import static com.hartwig.hmftools.cup.CuppaConfig.NOISE_ALLOCATIONS_DESC;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_COUNTS_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.REF_SNV_SAMPLE_POS_FREQ_FILE;
import static com.hartwig.hmftools.cup.CuppaConfig.configCategories;
import static com.hartwig.hmftools.cup.prep.PrepConfig.addPipelineDirectories;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.common.NoiseRefCache;
import com.hartwig.hmftools.cup.feature.RefFeatures;
import com.hartwig.hmftools.cup.rna.RefGeneExpression;
import com.hartwig.hmftools.cup.somatics.RefSomatics;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class RefDataConfig
{
    public final RefGenomeVersion RefGenVersion;
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
    public final String VirusDir;
    public final String IsofoxDir;
    public final String SomaticVariantsDir;

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
    public static final String SOMATIC_VARIANTS_DIR = "somatic_variants_dir";

    private static final String REF_FEATURE_OVERRIDE_FILE = "feature_override_file";

    public static final String GENDER_RATES = "gender_rates";
    public static final String GENDER_RATES_ADULT_DEFAULT = "ADULT_DEFAULT";

    private static final String WRITE_COHORT_FILES = "write_cohort_files";

    private static final String FILE_DELIM = ",";

    public RefDataConfig(final ConfigBuilder configBuilder)
    {
        Categories = configCategories(configBuilder);

        CUP_LOGGER.info("build ref data for classifiers: {}", Categories.isEmpty() ? ALL_CATEGORIES : Categories.toString());

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        SampleDataFile = configBuilder.getValue(REF_SAMPLE_DATA_FILE, "");
        CohortSampleTraitsFile = configBuilder.getValue(REF_COHORT_SAMPLE_TRAITS_FILE, "");
        CohortSigContribsFile = configBuilder.getValue(REF_COHORT_SIG_CONTRIBS_FILE, "");
        CohortSampleSvDataFile = configBuilder.getValue(REF_COHORT_SV_DATA_FILE, "");
        CohortFeaturesFile = configBuilder.getValue(REF_COHORT_FEATURES_FILE, "");

        LinxDir = configBuilder.getValue(LINX_DIR_CFG, "");
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, "");
        VirusDir = configBuilder.getValue(VIRUS_DIR_CFG, "");
        IsofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG, "");
        SomaticVariantsDir = configBuilder.getValue(SOMATIC_VARIANTS_DIR, "");

        GenPosMatrixFile = configBuilder.getValue(REF_SNV_SAMPLE_POS_FREQ_FILE, "");
        Snv96MatrixFile = configBuilder.getValue(REF_SNV_COUNTS_FILE, "");

        GeneExpMatrixFile = configBuilder.getValue(REF_GENE_EXP_DATA_FILE, "");
        AltSjMatrixFile = configBuilder.getValue(REF_ALT_SJ_DATA_FILE, "");

        FeatureOverrideFile = configBuilder.getValue(REF_FEATURE_OVERRIDE_FILE, "");

        DbAccess = createDatabaseAccess(configBuilder);

        OutputDir = parseOutputDir(configBuilder);

        NoiseAdjustments = new NoiseRefCache(OutputDir);
        NoiseAdjustments.loadNoiseAllocations(configBuilder.getValue(NOISE_ALLOCATIONS));

        WriteCohortFiles = configBuilder.hasFlag(WRITE_COHORT_FILES);
    }

    public static final List<String> parseFileSet(final String filenames)
    {
        return Arrays.stream(filenames.split(FILE_DELIM, -1)).collect(Collectors.toList());
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(CATEGORIES, false, "Categories to build ref data for");

        configBuilder.addPath(REF_SAMPLE_DATA_FILE, true, "Ref sample data file");

        // when combining datasets, these paths are ',' delimited so cannot be checked as a path
        // could add logic for this to ConfigBuilder but would be rarely used
        configBuilder.addConfigItem(REF_SNV_SAMPLE_POS_FREQ_FILE, false, "Ref SNV position frequency matrix data file");
        configBuilder.addConfigItem(REF_SNV_COUNTS_FILE, false, "Ref SNV trinucleotide matrix data file");
        configBuilder.addConfigItem(REF_COHORT_FEATURES_FILE, false, "Ref sample features data file");
        configBuilder.addConfigItem(REF_COHORT_SAMPLE_TRAITS_FILE, false, "Ref sample cohort traits file");
        configBuilder.addConfigItem(REF_COHORT_SIG_CONTRIBS_FILE, false, "Ref sample cohort signature contributions file");
        configBuilder.addConfigItem(REF_COHORT_SV_DATA_FILE, false, "Ref sample cohort SV file");

        addPipelineDirectories(configBuilder);

        configBuilder.addPath(REF_GENE_EXP_DATA_FILE, false, "Ref sample RNA gene expression cohort data file");
        configBuilder.addPath(REF_ALT_SJ_DATA_FILE, false, "Ref sample RNA alternate splice junction cohort data file");
        configBuilder.addPath(SOMATIC_VARIANTS_DIR, false, "Directory with flat file (converted) somatic variant files");

        configBuilder.addPath(REF_FEATURE_OVERRIDE_FILE, false, "Ref feature override data file");

        configBuilder.addConfigItem(NOISE_ALLOCATIONS, false, NOISE_ALLOCATIONS_DESC);
        configBuilder.addConfigItem(GENDER_RATES, false, "Gender-rate overrides - format CancerType;MalePerc;FemalePerc, etc");
        configBuilder.addFlag(WRITE_COHORT_FILES, "Re-write ref data as cohort files");

        addDatabaseCmdLineArgs(configBuilder, false);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);

        RefGeneExpression.registerConfig(configBuilder);
        RefSomatics.registerConfig(configBuilder);
    }
}
