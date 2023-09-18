package com.hartwig.hmftools.sieve;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_REGIONS;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificRegions;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class MapDropperConfig
{
    public static final Logger MD_LOGGER = LogManager.getLogger(MapDropperConfig.class);
    private static final String BAM_FILE = "bam_file";
    private static final String OUTPUT_BAM_FILE = "output_bam_file";
    public final String BamFile;
    public final String RefGenomeFile;
    public final String OutputBamFile;
    public final RefGenomeVersion RefGenVersion;
    public final List<ChrBaseRegion> SpecificRegions;
    public final int Threads;

    public MapDropperConfig(@NotNull final ConfigBuilder configBuilder)
    {
        BamFile = configBuilder.getValue(BAM_FILE);
        OutputBamFile = configBuilder.getValue(OUTPUT_BAM_FILE);
        RefGenomeFile = configBuilder.getValue(REF_GENOME);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
        Threads = parseThreads(configBuilder);

        SpecificRegions = Lists.newArrayList();
        try
        {
            SpecificRegions.addAll(loadSpecificRegions(configBuilder.getValue(SPECIFIC_REGIONS)));
        }
        catch(ParseException e)
        {
            MD_LOGGER.error("failed to load specific regions");
        }
    }

    public static void addConfig(@NotNull final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(BAM_FILE, true, "BAM file to drop bad maps from");
        configBuilder.addConfigItem(OUTPUT_BAM_FILE, true, "Output BAM file");
        addRefGenomeConfig(configBuilder, true);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);
    }
}
