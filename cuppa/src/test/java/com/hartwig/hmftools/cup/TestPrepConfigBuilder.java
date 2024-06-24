package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;

import java.util.Arrays;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class TestPrepConfigBuilder
{
    public static final String TEST_SAMPLE_ID_FILE = Resources.getResource("pipeline_output/sample_ids.csv").getPath();
    public static final List<String> TEST_SAMPLE_IDS = Arrays.asList("COLO829v003T");
    public static final RefGenomeVersion TEST_REF_GENOME_VERSION = RefGenomeVersion.V37;
    public static final List<CategoryType> TEST_CATEGORIES = CategoryType.getDnaCategories();
    public static final String TEST_OUTPUT_DIR = "";
    public static final String TEST_OUTPUT_ID = null;
    public static final Integer TEST_THREADS = 1;
    public static final boolean TEST_WRITE_BY_CATEGORY = true;
    public static final String TEST_ALT_SPLICE_JUNCTION_SITES = Resources.getResource("alt_sj.selected_loci.minimal.tsv").getPath();
    public static final String TEST_SAMPLE_DATA_DIR = Resources.getResource("pipeline_output/").getPath();
    public static final String TEST_SOMATIC_VARIANTS_DIR = Resources.getResource("liftover_output/").getPath();

    private List<String> SampleIds = TEST_SAMPLE_IDS;
    public RefGenomeVersion RefGenVersion = TEST_REF_GENOME_VERSION;
    public List<CategoryType> Categories = TEST_CATEGORIES;
    private String OutputDir = TEST_OUTPUT_DIR;
    private String OutputId = TEST_OUTPUT_ID; // for multi-sample mode
    private int Threads = TEST_THREADS;
    private boolean WriteByCategory = TEST_WRITE_BY_CATEGORY;

    private String SampleDataDir = "";
    private String LinxDir = SampleDataDir; // pipeline directories, accepting wildcards
    private String PurpleDir = SampleDataDir;
    private String VirusDir = SampleDataDir;
    private String IsofoxDir = SampleDataDir;
    private String SomaticVariantsDir = TEST_SOMATIC_VARIANTS_DIR;

    private String AltSpliceJunctionSites = TEST_ALT_SPLICE_JUNCTION_SITES;
    private int ProgressInterval = 100;


    public TestPrepConfigBuilder sampleIds(List<String> sampleIds)
    {
        SampleIds = sampleIds;
        return this;
    }

    public TestPrepConfigBuilder refGenomeVersion(RefGenomeVersion refGenomeVersion)
    {
        RefGenVersion = refGenomeVersion;
        return this;
    }

    public TestPrepConfigBuilder refGenomeVersion(String refGenomeVersion)
    {
        RefGenVersion = RefGenomeVersion.valueOf(refGenomeVersion);
        return this;
    }

    public TestPrepConfigBuilder categories(List<CategoryType> categories)
    {
        Categories = categories;
        return this;
    }

    public TestPrepConfigBuilder outputDir(String outputDir)
    {
        OutputDir = outputDir;
        return this;
    }

    public TestPrepConfigBuilder threads(int threads)
    {
        Threads = threads;
        return this;
    }

    public TestPrepConfigBuilder outputId(String outputId)
    {
        OutputId = outputId;
        return this;
    }

    public TestPrepConfigBuilder writeByCategory(boolean writeByCategory)
    {
        WriteByCategory = writeByCategory;
        return this;
    }

    public TestPrepConfigBuilder sampleDataDir(String sampleDataDir)
    {
        SampleDataDir = sampleDataDir;
        PurpleDir = sampleDataDir;
        LinxDir = sampleDataDir;
        VirusDir = sampleDataDir;
        IsofoxDir = sampleDataDir;
        SomaticVariantsDir = sampleDataDir;
        return this;
    }

    public TestPrepConfigBuilder linxDir(String linxDir)
    {
        LinxDir = linxDir;
        return this;
    }

    public TestPrepConfigBuilder purpleDir(String purpleDir)
    {
        PurpleDir = purpleDir;
        return this;
    }

    public TestPrepConfigBuilder virusDir(String virusDir)
    {
        VirusDir = virusDir;
        return this;
    }

    public TestPrepConfigBuilder isofoxDir(String isofoxDir)
    {
        IsofoxDir = isofoxDir;
        return this;
    }

    public TestPrepConfigBuilder somaticVariantsDir(String somaticVariantsDir)
    {
        SomaticVariantsDir = somaticVariantsDir;
        return this;
    }

    public TestPrepConfigBuilder altSpliceJunctionSites(String altSpliceJunctionSites)
    {
        AltSpliceJunctionSites = altSpliceJunctionSites;
        return this;
    }

    public TestPrepConfigBuilder progressInterval(int progressInterval)
    {
        ProgressInterval = progressInterval;
        return this;
    }

    public PrepConfig build()
    {
        return new PrepConfig(
                SampleIds,
                Categories,
                RefGenVersion,
                OutputDir,
                OutputId,
                Threads,
                WriteByCategory,
                SampleDataDir,
                LinxDir,
                PurpleDir,
                VirusDir,
                IsofoxDir,
                SomaticVariantsDir,
                AltSpliceJunctionSites,
                ProgressInterval
        );
    };

    public static PrepConfig fromArgs(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PrepConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        return new PrepConfig(configBuilder);
    }
}
