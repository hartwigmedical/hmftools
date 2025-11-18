package com.hartwig.hmftools.bamtools.unmappableregions;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class UnmappableRegionsConfig
{
    public final String RefGenome;
    public final boolean IncludeAllChromosomes;

    public final String OutputBedFile;

    private static final String INCLUDE_ALL_CHROMOSOMES = "all_chromosomes";
    private static final String INCLUDE_ALL_CHROMOSOMES_DESC = "Include all chromosomes (including non-human, e.g. MT)";

    private static final String OUTPUT_BED = "output_bed";
    private static final String OUTPUT_BED_DESC = "Output BED file";

    public UnmappableRegionsConfig(final ConfigBuilder configBuilder)
    {
        RefGenome = configBuilder.getValue(REF_GENOME);
        IncludeAllChromosomes = configBuilder.hasFlag(INCLUDE_ALL_CHROMOSOMES);
        OutputBedFile = configBuilder.getValue(OUTPUT_BED);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        RefGenomeSource.addRefGenomeFile(configBuilder, true);
        configBuilder.addConfigItem(OUTPUT_BED, true, OUTPUT_BED_DESC);
        configBuilder.addFlag(INCLUDE_ALL_CHROMOSOMES, INCLUDE_ALL_CHROMOSOMES_DESC);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
