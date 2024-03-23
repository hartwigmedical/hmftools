package com.hartwig.hmftools.cobalt.metrics;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.from;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MetricsConfig
{
    public final String SampleId;
    public final String BamFile;
    public final String RefGenomeFile;
    public final RefGenomeVersion RefGenVersion;
    public final RefGenomeInterface RefGenome;

    public final int PartitionSize;
    public final Map<String,List<BaseRegion>> TargetRegions;
    public final boolean OnlyTargetRegions;
    public final String OutputDir;
    public final String OutputId;

    public final int FragmentLengthUnits;
    public final double GcPercentUnits;
    public final boolean CaptureRegionCounts;

    public final int Threads;

    public final SpecificRegions SpecificChrRegions;

    private static final String BAM_FILE = "bam_file";
    private static final String ONLY_TARGET = "only_target";
    private static final String FRAG_LENGTH_UNITS = "frag_length_units";
    private static final String GC_PERCENT_UNITS = "gc_percent_units";
    private static final String CAPTURE_REGION_COUNTS = "capture_region_counts";

    public static final int TARGET_REGION_PROXIMITY = 100;

    private static final int DEFAULT_PARTITION_SIZE = 1_000_000;
    private static final double DEFAULT_GC_BUCKET = 0.01;
    private static final int DEFAULT_FRAG_LENGTH_BUCKET = 10;

    public MetricsConfig(final ConfigBuilder configBuilder)
    {
        SampleId =  configBuilder.getValue(SAMPLE);
        BamFile =  configBuilder.getValue(BAM_FILE);
        RefGenomeFile =  configBuilder.getValue(REF_GENOME);

        RefGenome = RefGenomeSource.loadRefGenome(RefGenomeFile);

        if(configBuilder.hasValue(OUTPUT_DIR))
        {
            OutputDir = parseOutputDir(configBuilder);
        }
        else
        {
            OutputDir = checkAddDirSeparator(Paths.get(BamFile).getParent().toString());
        }

        OutputId =  configBuilder.getValue(OUTPUT_ID);

        RefGenVersion = configBuilder.hasValue(REF_GENOME_VERSION) ? from(configBuilder) : deriveRefGenomeVersion(BamFile);
        CB_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        CB_LOGGER.info("output({})", OutputDir);

        PartitionSize = DEFAULT_PARTITION_SIZE;

        TargetRegions = loadChrBaseRegions(configBuilder.getValue(TARGET_REGIONS_BED));
        OnlyTargetRegions = !TargetRegions.isEmpty() && configBuilder.hasFlag(ONLY_TARGET);

        FragmentLengthUnits = configBuilder.getInteger(FRAG_LENGTH_UNITS);
        GcPercentUnits = configBuilder.getDecimal(GC_PERCENT_UNITS);
        CaptureRegionCounts = configBuilder.hasFlag(CAPTURE_REGION_COUNTS);

        SpecificChrRegions = SpecificRegions.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        addRefGenomeVersion(configBuilder);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, true, "BAM file path");

        configBuilder.addPath(TARGET_REGIONS_BED, true, TARGET_REGIONS_BED_DESC);

        configBuilder.addInteger(FRAG_LENGTH_UNITS, "Fragment length bucket size", DEFAULT_FRAG_LENGTH_BUCKET);
        configBuilder.addDecimal(GC_PERCENT_UNITS, "GC percent bucket size", DEFAULT_GC_BUCKET);

        configBuilder.addFlag(ONLY_TARGET, "Only capture metrics within the specific regions file");
        configBuilder.addFlag(CAPTURE_REGION_COUNTS, "Capture metrics for each regions");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
