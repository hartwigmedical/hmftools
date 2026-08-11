package com.hartwig.hmftools.bamtools.synthetic;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BAM_FILE_DESC;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.REGIONS_FILE;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.DEFAULT_CN_BACKBONE_DISTANCE;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.DEFAULT_MAX_SMALL_VARIANTS;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.DEFAULT_MAX_SVS;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.bamops.BamToolName.BAMTOOL_PATH;
import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.loadDriverGenes;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;

import java.util.List;

import com.hartwig.hmftools.bamtools.slice.SliceParams;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SyntheticConfig
{
    public final String SampleId;
    public final String ReferenceId;
    public final String BamFile;
    public final String OutputPrefix;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;

    public final String RegionsFilename;
    public final String PurpleDir;
    public final String LinxDir;
    public final String LinxGermlineDir;
    public final List<DriverGene> DriverGenes;
    public final String EnsemblDataDir;
    public final String GermlineHetSites;
    public final String OutputDir;

    public final boolean WriteBam;
    public final boolean WriteRegionData;

    public final SliceParams Params;

    public final int Threads;
    public final int MaxSmallVariants;
    public final int MaxSVs;
    public final int CopyNumberBackboneDistance;
    public final String BamToolPath;

    private static final String OUTPUT_PREFIX = "output_prefix";
    private static final String WRITE_BAM = "write_bam";
    private static final String WRITE_REGION_DATA = "write_region_data";

    private static final String MAX_SMALL_VARIANTS = "max_small_variants";
    private static final String MAX_SVS = "max_svs";
    private static final String CN_BACKBONE_DISTANCE = "cn_backbone_distance";
    private static final String GERMLINE_HET_SITES = "germline_het_sites";

    public SyntheticConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        ReferenceId = configBuilder.getValue(REFERENCE);
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : pathFromFile(BamFile);

        if(configBuilder.hasValue(OUTPUT_PREFIX))
        {
            OutputPrefix = configBuilder.getValue(OUTPUT_PREFIX);
        }
        else
        {
            String filename = filenamePart(BamFile);
            int extIndex = filename.lastIndexOf(".");
            String bamFileName = filename.substring(0, extIndex);

            OutputPrefix = bamFileName + ".slice";
        }

        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        LinxDir = configBuilder.getValue(LINX_DIR_CFG);
        LinxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR_CFG, LinxDir);

        DriverGenes = loadDriverGenes(configBuilder);
        EnsemblDataDir = configBuilder.getValue(EnsemblDataCache.ENSEMBL_DATA_DIR);
        GermlineHetSites = configBuilder.getValue(GERMLINE_HET_SITES);

        BamToolPath = configBuilder.getValue(BAMTOOL_PATH);

        WriteRegionData = configBuilder.hasFlag(WRITE_REGION_DATA);
        WriteBam = configBuilder.hasFlag(WRITE_BAM) || configBuilder.hasValue(BAMTOOL_PATH);

        CopyNumberBackboneDistance = configBuilder.getInteger(CN_BACKBONE_DISTANCE);
        MaxSmallVariants = configBuilder.getInteger(MAX_SMALL_VARIANTS);
        MaxSVs = configBuilder.getInteger(MAX_SVS);

        if(BamFile == null || OutputDir == null || RefGenomeFile == null)
        {
            BT_LOGGER.error("missing config: bam({}) refGenome({}) outputDir({})",
                    BamFile != null, RefGenomeFile != null, OutputDir != null);
            System.exit(1);
        }

        RefGenVersion = deriveRefGenomeVersion(BamFile);
        BT_LOGGER.info("input bam({}) refGenome({})", BamFile, RefGenVersion);
        BT_LOGGER.info("output({}) and file prefix({})", OutputDir, OutputPrefix);

        RegionsFilename = configBuilder.getValue(REGIONS_FILE);

        Params = new SliceParams(configBuilder);

        Threads = parseThreads(configBuilder);
    }

    public boolean isPanelMode() { return RegionsFilename != null; }
    public boolean isWgsMode() { return !isPanelMode(); }

    /*
    public String formFilename(final String fileExtension)
    {
        String outputFile = OutputDir + SampleId + "." + OutputPrefix + fileExtension;
        return outputFile;
    }
    */

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);
        configBuilder.addPath(BAM_FILE, true, BAM_FILE_DESC);

        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_CFG);
        configBuilder.addPath(REGIONS_FILE, false, "Region or panel definition file");
        configBuilder.addPath(GERMLINE_HET_SITES, false, "Germline heterozygous sites");

        DriverGenePanelConfig.addGenePanelOption(configBuilder, false);
        EnsemblDataCache.addEnsemblDir(configBuilder);

        addRefGenomeFile(configBuilder, true);;
        configBuilder.addConfigItem(OUTPUT_PREFIX, true,"File prefix for BAM and region info TSV");

        configBuilder.addInteger(CN_BACKBONE_DISTANCE, "Max structural variants", DEFAULT_CN_BACKBONE_DISTANCE);
        configBuilder.addInteger(MAX_SMALL_VARIANTS, "Max structural variants", DEFAULT_MAX_SMALL_VARIANTS);
        configBuilder.addInteger(MAX_SVS, "Max structural variants", DEFAULT_MAX_SVS);

        configBuilder.addFlag(WRITE_BAM, "Write BAM file for sliced region");
        configBuilder.addFlag(WRITE_REGION_DATA, "Write slice region data");
        SliceParams.addConfig(configBuilder);

        BamToolName.addConfig(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }
}
