package com.hartwig.hmftools.redux.ms_sites;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE;
import static com.hartwig.hmftools.common.genome.gc.GCProfileFactory.GC_PROFILE_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MsFinderConfig
{
    public final String RefGenomeFile;
    public final String OutputDir;

    public final String GcProfilePath;
    public final String EnsemblCacheDir;

    public final int Threads;
    public final SpecificRegions SpecificChrRegions;

    // constants & config
    protected static final int DOWNSAMPLE_FACTOR = 25_000_000;
    protected static final int BED_REGION_EXPANSION = 950;
    protected static final int MIN_TARGET_SITE_COUNT = 15_000;
    protected static final int CHUNK_SIZE = 100_000;
    protected static final double MIN_MAPPABILITY = 0.7;

    public MsFinderConfig(final ConfigBuilder configBuilder)
    {
        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        GcProfilePath = configBuilder.getValue(GC_PROFILE);
        EnsemblCacheDir = configBuilder.getValue(ENSEMBL_DATA_DIR);
        Threads = parseThreads(configBuilder);
        OutputDir = parseOutputDir(configBuilder);

        SpecificChrRegions = SpecificRegions.from(configBuilder);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(GC_PROFILE, true, GC_PROFILE_DESC);
        configBuilder.addPath("bed_file", false, "exta bed file input");

        addOutputDir(configBuilder);

        addEnsemblDir(configBuilder, true);
        addSpecificChromosomesRegionsConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);

        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
