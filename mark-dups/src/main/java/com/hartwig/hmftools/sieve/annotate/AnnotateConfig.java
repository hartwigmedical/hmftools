package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class AnnotateConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(AnnotateConfig.class);

    public final String BamFile;
    public final String BedFile;
    public final String BlacklistRegionRepeatMaskerFile;
    public final String RefGenome;
    public final String OutputFile;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    private static final String BAM_FILE = "bam_file";
    private static final String BED_FILE = "bed_file";
    private static final String BLACKLIST_REGION_REPEAT_MASKER_FILE = "blacklist_region_repeat_makser_file";

    private static final String OUTPUT_FILE = "output_file";

    public AnnotateConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        BedFile = configBuilder.getValue(BED_FILE);
        BlacklistRegionRepeatMaskerFile = configBuilder.getValue(BLACKLIST_REGION_REPEAT_MASKER_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file to annotate");
        configBuilder.addPath(BED_FILE, true, "BED file containing bad regions");
        configBuilder.addConfigItem(BLACKLIST_REGION_REPEAT_MASKER_FILE, true, "CSV file containing the blacklisted regions and their repeat masks");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
    }
}
