package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.perf.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.prep.SampleIdsLoader.COL_RNA_SAMPLE_ID;
import static com.hartwig.hmftools.cup.prep.SampleIdsLoader.COL_SAMPLE_ID;
import static com.hartwig.hmftools.cup.prep.SampleIdsLoader.NO_RNA_SAMPLE_ID;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.SOMATIC_VARIANTS_DIR_CFG;
import static com.hartwig.hmftools.cup.somatics.SomaticVariant.SOMATIC_VARIANTS_DIR_DESC;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.cup.somatics.SomaticVariant;

public class PrepConfig
{
    public final List<String> SampleIds;
    public final List<String> RnaSampleIds;

    // pipeline directories, accepting wildcards
    public final String SampleDataDir;
    public final String LinxDir;
    public final String PurpleDir;
    public final String VirusDir;
    public final String IsofoxDir;
    public final String SomaticVariantsDir;

    public final List<CategoryType> Categories;
    public final RefGenomeVersion RefGenVersion;
    public final String AltSpliceJunctionSites;

    public final String OutputDir;
    public final String OutputId; // for multi-sample mode

    public final boolean WriteByCategory;
    public final int Threads;

    public static final String CATEGORIES = "categories";
    public static final String CATEGORIES_DESC = "Categories to build ref data for";
    public static final String ALL_CATEGORIES = "ALL";
    public static final String DNA_CATEGORIES = "DNA";
    public static final String RNA_CATEGORIES = "RNA";
    public static final String CATEGORY_DELIM = ";";

    public static final String REF_ALT_SJ_SITES = "ref_alt_sj_sites";
    public static final String REF_ALT_SJ_SITES_DESC = "RNA required alternative splice junction sites";

    public static final String WRITE_FILE_BY_CATEGORY = "write_by_category";
    public static final String WRITE_FILE_BY_CATEGORY_DESC = "Cohort mode - write files by category";

    public static final String THREADS_DESC = "Number of threads to use in multi sample mode";

    public static final String RNA_SAMPLE_ID_DESC = String.format("RNA sample ID. Note: Output file prefix uses value from `-%s`", SAMPLE);

    public static final String SAMPLE_ID_FILE_DESC =
            String.format("TSV file with columns %s and/or %s. ", COL_SAMPLE_ID, COL_RNA_SAMPLE_ID) +
            String.format("%s values can be blank (%s value used instead) or %s (skips extracting RNA features)", COL_RNA_SAMPLE_ID, COL_SAMPLE_ID, NO_RNA_SAMPLE_ID);

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        SampleIdsLoader sampleIdsLoader = new SampleIdsLoader().loadFromConfig(configBuilder);
        SampleIds = sampleIdsLoader.sampleIds();
        RnaSampleIds = sampleIdsLoader.rnaSampleIds();

        Categories = parseCategories(configBuilder);

        SampleDataDir = configBuilder.getValue(SAMPLE_DATA_DIR_CFG, "");
        LinxDir = configBuilder.getValue(LINX_DIR_CFG, SampleDataDir);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir);
        VirusDir = configBuilder.getValue(VIRUS_DIR_CFG, SampleDataDir);
        IsofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG, SampleDataDir);
        SomaticVariantsDir = configBuilder.getValue(SOMATIC_VARIANTS_DIR_CFG, SampleDataDir);

        RefGenVersion = RefGenomeVersion.from(configBuilder);
        AltSpliceJunctionSites = configBuilder.getValue(REF_ALT_SJ_SITES);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        WriteByCategory = configBuilder.hasFlag(WRITE_FILE_BY_CATEGORY);

        Threads = TaskExecutor.parseThreads(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, RNA_SAMPLE_ID_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_CFG);
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
        configBuilder.addPath(SOMATIC_VARIANTS_DIR_CFG, false, SOMATIC_VARIANTS_DIR_DESC);

        configBuilder.addConfigItem(CATEGORIES, false, CATEGORIES_DESC);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addPath(REF_ALT_SJ_SITES, false, REF_ALT_SJ_SITES_DESC);

        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        configBuilder.addConfigItem(OUTPUT_ID, false, OUTPUT_ID_DESC);

        configBuilder.addFlag(WRITE_FILE_BY_CATEGORY, WRITE_FILE_BY_CATEGORY_DESC);
        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");

        ConfigUtils.addLoggingOptions(configBuilder);
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public boolean isSingleSample() { return SampleIds.size() == 1; }

    // Generate input file paths by sample id
    public String getLinxDataDir(final String sampleId) { return ConfigUtils.convertWildcardSamplePath(LinxDir, sampleId); }
    public String getPurpleDataDir(final String sampleId) { return ConfigUtils.convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getVirusDataDir(final String sampleId) { return ConfigUtils.convertWildcardSamplePath(VirusDir, sampleId); }
    public String getIsofoxDataDir(final String sampleId) { return ConfigUtils.convertWildcardSamplePath(IsofoxDir, sampleId); }
    public String getSomaticVariantsDir(final String sampleId) { return ConfigUtils.convertWildcardSamplePath(SomaticVariantsDir, sampleId); }

    public String purpleSomaticVcfFile(final String sampleId) { return PurpleCommon.purpleSomaticVcfFile(getPurpleDataDir(sampleId), sampleId); }
    public String somaticVariantsGenericFile(final String sampleId) { return SomaticVariant.generateFilename(getSomaticVariantsDir(sampleId), sampleId); }
    public String purpleSvFile(final String sampleId) { return PurpleCommon.purpleSomaticSvFile(getPurpleDataDir(sampleId), sampleId); }
    public String purplePurityFile(final String sampleId) { return PurpleCommon.purplePurityFile(getPurpleDataDir(sampleId), sampleId); }
    public String purpleDriverCatalogFile(final String sampleId) { return DriverCatalogFile.generateSomaticFilename(getPurpleDataDir(sampleId), sampleId); }
    public String linxDriverCatalogFile(final String sampleId) { return LinxDriver.generateCatalogFilename(getLinxDataDir(sampleId), sampleId, true); }
    public String linxClusterFile(final String sampleId) { return LinxCluster.generateFilename(getLinxDataDir(sampleId), sampleId, false); }
    public String linxFusionFile(final String sampleId) { return LinxFusion.generateFilename(getLinxDataDir(sampleId), sampleId); }
    public String viralAnnotationFile(final String sampleId) { return AnnotatedVirusFile.generateFileName(getVirusDataDir(sampleId), sampleId); }
    public String geneExpressionFile(final String sampleId) { return GeneExpressionFile.generateFilename(getIsofoxDataDir(sampleId), sampleId); }
    public String altSpliceJunctionFile(final String sampleId) { return AltSpliceJunctionFile.generateFilename(getIsofoxDataDir(sampleId), sampleId); }

    private static List<CategoryType> parseCategories(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(CATEGORIES))
            return CategoryType.getDnaCategories(); // Default to DNA

        String configCategories = configBuilder.getValue(CATEGORIES).toUpperCase();

        switch(configCategories)
        {
            case ALL_CATEGORIES ->
            {
                return CategoryType.getAllCategories();
            }
            case DNA_CATEGORIES ->
            {
                return CategoryType.getDnaCategories();
            }
            case RNA_CATEGORIES ->
            {
                return CategoryType.getRnaCategories();
            }
        }

        final String[] categoryStrings = configCategories.split(CATEGORY_DELIM);
        return Arrays.stream(categoryStrings).map(CategoryType::valueOf).collect(Collectors.toList());
    }

    @VisibleForTesting
    public PrepConfig(
            final List<String> sampleIds,
            final List<String> rnaSampleIds,
            final List<CategoryType> categories,
            final RefGenomeVersion refGenVersion,
            final String outputDir,
            final String outputId,
            final int threads,
            final boolean writeByCategory,
            final String sampleDataDir,
            final String linxDir,
            final String purpleDir,
            final String virusDir,
            final String isofoxDir,
            final String somaticVariantsDir,
            final String altSpliceJunctionSites
    )
    {
        SampleIds = sampleIds;
        RnaSampleIds = rnaSampleIds;
        Categories = categories;
        RefGenVersion = refGenVersion;
        OutputDir = outputDir;
        OutputId = outputId;
        Threads = threads;
        WriteByCategory = writeByCategory;
        SampleDataDir = sampleDataDir;
        LinxDir = linxDir;
        PurpleDir = purpleDir;
        VirusDir = virusDir;
        IsofoxDir = isofoxDir;
        SomaticVariantsDir = somaticVariantsDir;
        AltSpliceJunctionSites = altSpliceJunctionSites;
    }
}
