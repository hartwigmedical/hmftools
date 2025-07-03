package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class PanelConfig
{
    public final String AmberSitesFile;
    public final String GeneTranscriptFile;
    public final String CustomRegionFile;

    public final String OutputPrefix;
    public final String OutputDir;

    public final RefGenomeInterface RefGenome;
    public final RefGenomeVersion RefGenVersion;

    public final int Threads;
    public final String EnsemblDir;

    private static final String AMBER_SITES_FILE = "amber_sites_file";
    private static final String GENE_TRANSCRIPT_FILE = "gene_transcript_file";
    private static final String CUSTOM_REGION_FILE = "custom_region_file";

    private static final String OUTPUT_PREFIX = "output_prefix";

    public PanelConfig(final ConfigBuilder configBuilder)
    {
        AmberSitesFile = configBuilder.getValue(AMBER_SITES_FILE);
        GeneTranscriptFile = configBuilder.getValue(GENE_TRANSCRIPT_FILE);
        CustomRegionFile = configBuilder.getValue(CUSTOM_REGION_FILE);

        final String refGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFile);
        RefGenVersion = deriveRefGenomeVersion((RefGenomeSource)RefGenome);

        EnsemblDir = configBuilder.getValue(EnsemblDataCache.ENSEMBL_DATA_DIR);

        OutputPrefix = configBuilder.getValue(OUTPUT_PREFIX);
        OutputDir = parseOutputDir(configBuilder);
        Threads = parseThreads(configBuilder);
    }

    public String formOutputFilename(final String fileIdExtension)
    {
        return OutputDir + OutputPrefix + fileIdExtension;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(AMBER_SITES_FILE, false, "Amber het sites file");
        configBuilder.addPath(GENE_TRANSCRIPT_FILE, false, "Gene and transcript name file");
        configBuilder.addPath(CUSTOM_REGION_FILE, false, "Custom region file");

        EnsemblDataCache.addEnsemblDir(configBuilder, true);

        addRefGenomeFile(configBuilder, true);

        configBuilder.addConfigItem(OUTPUT_PREFIX, true, "prefix of output BED filename");
        addOutputDir(configBuilder);

        addThreadOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
