package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.configCategories;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;

public class PrepConfig
{
    public final RefGenomeVersion RefGenVersion;
    public final List<CategoryType> Categories;

    public final List<String> SampleIds;

    public final String OutputDir;
    public final String OutputId; // for multi-sample mode
    public final int Threads;
    public final boolean WriteByCategory;

    // pipeline directories, accepting wildcards
    public final String SampleDataDir;
    public final String LinxDir;
    public final String PurpleDir;
    public final String VirusDir;
    public final String IsofoxDir;
    public final String AltSpliceJunctionSites;

    private static final String WRITE_FILE_BY_CATEGORY = "write_by_category";
    private static final String REF_ALT_SJ_SITES = "ref_alt_sj_sites";

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        Categories = configCategories(configBuilder);

        SampleIds = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            SampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            SampleIds.addAll(loadSampleIdsFile(configBuilder));
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        SampleDataDir = configBuilder.getValue(SAMPLE_DATA_DIR_CFG, "");
        LinxDir = configBuilder.getValue(LINX_DIR_CFG, SampleDataDir);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, SampleDataDir);
        VirusDir = configBuilder.getValue(VIRUS_DIR_CFG, SampleDataDir);
        IsofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG, SampleDataDir);

        AltSpliceJunctionSites = configBuilder.getValue(REF_ALT_SJ_SITES);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
        Threads = parseThreads(configBuilder);

        WriteByCategory = SampleIds.size() > 1 && configBuilder.hasFlag(WRITE_FILE_BY_CATEGORY);
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public boolean isSingleSample() { return SampleIds.size() == 1; }

    // Add args to config builder
    public static void addPipelineDirectories(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, VIRUS_DIR_CFG);
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(CATEGORIES, false, "Categories to build ref data for");

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addPath(REF_ALT_SJ_SITES, false, "RNA required alternative splice junction sites");
        configBuilder.addFlag(WRITE_FILE_BY_CATEGORY, "Cohort mode - write files by category");
        configBuilder.addConfigItem(THREADS, false, "Number of threads to use in multi sample mode", "1");

        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);
        addPipelineDirectories(configBuilder);

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
    }

    // Generate input file paths by sample id
    public String getLinxDataDir(final String sampleId) { return convertWildcardSamplePath(LinxDir, sampleId); }
    public String getPurpleDataDir(final String sampleId) { return convertWildcardSamplePath(PurpleDir, sampleId); }
    public String getVirusDataDir(final String sampleId) { return convertWildcardSamplePath(VirusDir, sampleId); }
    public String getIsofoxDataDir(final String sampleId) { return convertWildcardSamplePath(IsofoxDir, sampleId); }

    public String purpleSomaticVcfFile(final String sampleId) { return PurpleCommon.purpleSomaticVcfFile(getPurpleDataDir(sampleId), sampleId); }
    public String purpleSvFile(final String sampleId) { return PurpleCommon.purpleSomaticSvFile(getPurpleDataDir(sampleId), sampleId); }
    public String purplePurityFile(final String sampleId) { return PurpleCommon.purplePurityFile(getPurpleDataDir(sampleId), sampleId); }
    public String purpleQcFile(final String sampleId) { return PurpleCommon.purpleQcFile(getPurpleDataDir(sampleId), sampleId); }
    public String purpleDriverCatalogFile(final String sampleId) { return DriverCatalogFile.generateSomaticFilename(getPurpleDataDir(sampleId), sampleId); }
    public String linxDriverCatalogFile(final String sampleId) { return LinxDriver.generateCatalogFilename(getLinxDataDir(sampleId), sampleId, true); }
    public String linxClusterFile(final String sampleId) { return LinxCluster.generateFilename(getLinxDataDir(sampleId), sampleId, false); }
    public String linxFusionFile(final String sampleId) { return LinxFusion.generateFilename(getLinxDataDir(sampleId), sampleId); }
    public String viralAnnotationFile(final String sampleId) { return AnnotatedVirusFile.generateFileName(getVirusDataDir(sampleId), sampleId); }
    public String geneExpressionFile(final String sampleId) { return GeneExpressionFile.generateFilename(getIsofoxDataDir(sampleId), sampleId); }
    public String altSpliceJunctionFile(final String sampleId) { return AltSpliceJunctionFile.generateFilename(getIsofoxDataDir(sampleId), sampleId); }

    @VisibleForTesting
    public PrepConfig(
            final List<String> sampleIds,
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
            final String altSpliceJunctionSites
    )
    {
        SampleIds = sampleIds;
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
        AltSpliceJunctionSites = altSpliceJunctionSites;
    }
}
