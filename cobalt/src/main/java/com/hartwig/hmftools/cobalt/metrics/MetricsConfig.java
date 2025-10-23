package com.hartwig.hmftools.cobalt.metrics;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.bam.BamUtils.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.from;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.ChromosomeData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MetricsConfig
{
    public final String SampleId;
    public final String BamFile;
    public final RefGenomeVersion RefGenVersion;
    public final int PartitionSize;
    public final int WindowSize;
    public final String OutputDir;
    public final int Threads;
    public final SpecificRegions SpecificChrRegions;

    public static final String BAM_FILE = "bam_file";
    public static final String WINDOW_SIZE = "window_size";
    private static final int DEFAULT_PARTITION_SIZE = 1_000_000;
    private static final int DEFAULT_WINDOW_SIZE = 1000;
    public static final int MIN_MAPPING_QUALITY = 50;

    public MetricsConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        BamFile = configBuilder.getValue(BAM_FILE);
        RefGenVersion = configBuilder.hasValue(REF_GENOME_VERSION) ? from(configBuilder) : deriveRefGenomeVersion(BamFile);
        OutputDir = parseOutputDir(configBuilder);
        PartitionSize = DEFAULT_PARTITION_SIZE;
        WindowSize = configBuilder.getInteger(WINDOW_SIZE);

        CB_LOGGER.info("refGenome({}), bam({})", RefGenVersion, BamFile);
        CB_LOGGER.info("output({})", OutputDir);
        SpecificChrRegions = SpecificRegions.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeVersion(configBuilder);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addPath(BAM_FILE, true, "BAM file path");
        addOutputDir(configBuilder, true);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        configBuilder.addInteger(WINDOW_SIZE, "window size", DEFAULT_WINDOW_SIZE);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }

    public List<ChromosomeData> chromosomes()
    {
        RefGenomeCoordinates coordinates = RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        List<Chromosome> chromosomes = new ArrayList<>();
        if(SpecificChrRegions.Chromosomes.isEmpty())
        {
            chromosomes.addAll(coordinates.Lengths.keySet().stream().sorted().toList());
        }
        else
        {
            chromosomes.addAll(SpecificChrRegions.Chromosomes.stream().map(HumanChromosome::fromString).sorted().toList());
        }

        return chromosomes.stream().map(contig ->
        {
            int length = coordinates.Lengths.get(contig);
            String name = RefGenVersion.versionedChromosome(contig);
            return new ChromosomeData(name, length);
        }).sorted((o1, o2) -> ContigComparator.INSTANCE.compare(o1.Name, o2.Name)).toList();
    }

    public ListMultimap<Chromosome, Partition> createPartitions()
    {
        return PartitionBuilder.partitionChromosomes(chromosomes(), PartitionSize, WindowSize);
    }
}
