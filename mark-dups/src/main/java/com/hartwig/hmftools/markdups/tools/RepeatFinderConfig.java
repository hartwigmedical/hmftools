package com.hartwig.hmftools.markdups.tools;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RepeatFinderConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(RepeatFinderConfig.class);

    private static final String OUTPUT_FILE = "output_file";
    private static final String BED_FILE = "bed_file";
    private static final String MIN_REPEAT_BASES = "min_repeat_bases";

    public final String OutputFile;
    public final String BedFile;
    public final int MinRepeatBases;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    private static final int DEFAULT_MIN_REPEAT_BASES = 30;

    public RepeatFinderConfig(final ConfigBuilder configBuilder)
    {
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        BedFile = (configBuilder.hasValue(BED_FILE)) ? configBuilder.getValue(BED_FILE) : null;
        MinRepeatBases = configBuilder.getInteger(MIN_REPEAT_BASES);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);

    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addPath(BED_FILE, false, "Annotate repeats with overlapping regions from the BED file. Note that overlapping regions in the BED file are merged.");
        configBuilder.addInteger(MIN_REPEAT_BASES, "Minimum length of a repeat to output (must be at least four)", DEFAULT_MIN_REPEAT_BASES);
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
    }
}
