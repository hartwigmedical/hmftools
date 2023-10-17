package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
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

    private static final String BAM_FILE = "bam_file";
    private static final String HIGH_DEPTH_FILE = "high_depth_file";
    private static final String MASKED_REGION_FILE = "masked_region_file";
    private static final String OUTPUT_FILE = "output_file";
    private static final String KEEP_DUPLICATES = "keep_duplicates";

    public final String BamFile;
    public final String HighDepthFile;
    public final String RefGenome;
    public final String OutputFile;
    public final String MaskedRegionFile;
    public final boolean KeepDuplicates;
    public final RefGenomeVersion RefGenVersion;
    public final int Threads;

    public AnnotateConfig(final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        HighDepthFile = configBuilder.getValue(HIGH_DEPTH_FILE);
        MaskedRegionFile = configBuilder.getValue(MASKED_REGION_FILE);
        OutputFile = configBuilder.getValue(OUTPUT_FILE);
        KeepDuplicates = configBuilder.hasFlag(KEEP_DUPLICATES);
        RefGenome = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file");
        configBuilder.addPath(HIGH_DEPTH_FILE, true, "high depth file");
        configBuilder.addPath(MASKED_REGION_FILE, true, "A file describing the masked regions in the reference genome");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file");
        configBuilder.addFlag(KEEP_DUPLICATES, "Include duplicates the in the counts");
        addEnsemblDir(configBuilder, true);
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
    }
}
