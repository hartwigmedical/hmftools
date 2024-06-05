package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.bam.BamUtils.addValidationStringencyOption;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.SpecificRegions.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.SpecificRegions.loadSpecificRegions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.commons.cli.ParseException;

import htsjdk.samtools.ValidationStringency;

public class JitterAnalyserConfig
{
    public final String SampleId;
    public final String BamPath;
    public final RefGenomeVersion RefGenVersion;
    public final String RefGenomeFile;

    public final String RefGenomeMsiFile;
    public final String OutputDir;

    public final int MinMappingQuality;

    public final int MaxSitesPerType;

    public final int PartitionSize;
    public final ValidationStringency BamStringency;
    public final int Threads;
    public final boolean WritePlots;

    public final List<ChrBaseRegion> SpecificRegions;

    public static final String JITTER_MSI_SITES_FILE = "ref_genome_msi_file";
    public static final String JITTER_MSI_SITES_FILE_DESC = "Path to ref genome MSI sitesΩ tsv";

    public static final String JITTER_MAX_SITES_PER_TYPE = "max_sites_per_type";
    public static final String JITTER_MAX_SITES_PER_TYPE_DESC = "Max number of sites per microsatellite unit / length type";

    private static final String MIN_MAP_QUALITY = "min_map_quality";
    private static final String PARTITION_SIZE = "partition_size";

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 50;
    public static final int DEFAULT_NUM_SITES_PER_TYPE = 5_000;
    public static final int DEFAULT_PARTITION_SIZE = 1_000_000;

    public JitterAnalyserConfig(final ConfigBuilder configBuilder) throws ParseException
    {
        SampleId = configBuilder.getValue(SAMPLE);
        BamPath = configBuilder.getValue("bam");
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenomeMsiFile = configBuilder.getValue(JITTER_MSI_SITES_FILE);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);
        MinMappingQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(JITTER_MAX_SITES_PER_TYPE);
        PartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        BamStringency = BamUtils.validationStringency(configBuilder);
        WritePlots = true;
        SpecificRegions = loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS));
    }

    public JitterAnalyserConfig(
            final String sampleId, final RefGenomeVersion refGenVersion, final String refGenomeFile,
            final String refGenomeMsiFile, final String outputDir, int minMappingQuality, int maxSitesPerType, int threads, boolean writePlots)
    {
        SampleId = sampleId;
        BamPath = null;
        RefGenVersion = refGenVersion;
        RefGenomeFile = refGenomeFile;
        RefGenomeMsiFile = refGenomeMsiFile;
        OutputDir = outputDir;
        MinMappingQuality = minMappingQuality;
        MaxSitesPerType = maxSitesPerType;
        PartitionSize = DEFAULT_PARTITION_SIZE;
        BamStringency = ValidationStringency.STRICT;
        Threads = threads;
        WritePlots = writePlots;
        SpecificRegions = null;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, "sample id");
        configBuilder.addPath("bam", true, "path to bam file");

        addRefGenomeVersion(configBuilder);
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC + ", required when using CRAM files");

        configBuilder.addPath(JITTER_MSI_SITES_FILE, true, JITTER_MSI_SITES_FILE_DESC);

        addOutputDir(configBuilder);

        configBuilder.addInteger(
                MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used", DEFAULT_MIN_MAPPING_QUALITY);

        configBuilder.addInteger(JITTER_MAX_SITES_PER_TYPE, JITTER_MAX_SITES_PER_TYPE_DESC, DEFAULT_NUM_SITES_PER_TYPE);

        configBuilder.addInteger(PARTITION_SIZE, "size of the partitions/jobs processed by worker threads", DEFAULT_PARTITION_SIZE);

        addThreadOptions(configBuilder);
        addValidationStringencyOption(configBuilder);
        addLoggingOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public boolean isValid()
    {
        checkCreateOutputDir(OutputDir);
        return true;
    }
}
