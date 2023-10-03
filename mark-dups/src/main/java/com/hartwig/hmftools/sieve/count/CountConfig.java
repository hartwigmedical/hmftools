package com.hartwig.hmftools.sieve.count;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CountConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(CountConfig.class);
    private static final String BAM_FILE = "bam_file";
    private static final String BUCKET_SIZE = "bucket_size";
    private static final String CONTIG_COUNT_OUTPUT_FILE = "contig_count_output_file";
    private static final String BUCKET_COUNT_OUTPUT_FILE = "bucket_count_output_file";

    public final String BamFile;
    public final int BucketSize;
    public final String RefGenome;
    public final String ContigCountOutputFile;
    public final String BucketCountOutputFile;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    public CountConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        BucketSize = configBuilder.getInteger(BUCKET_SIZE);
        ContigCountOutputFile = configBuilder.getValue(CONTIG_COUNT_OUTPUT_FILE);
        BucketCountOutputFile = configBuilder.getValue(BUCKET_COUNT_OUTPUT_FILE);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file");
        configBuilder.addInteger(BUCKET_SIZE, "The size of the buckets in units of bp.", 100000);
        configBuilder.addConfigItem(CONTIG_COUNT_OUTPUT_FILE, true, "Output TSV for read counts per contig");
        configBuilder.addConfigItem(BUCKET_COUNT_OUTPUT_FILE, true, "Output TSV for read counts per bp bucket");
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
    }
}
