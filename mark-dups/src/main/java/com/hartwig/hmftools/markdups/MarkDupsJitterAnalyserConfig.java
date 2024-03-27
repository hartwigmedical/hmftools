package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConfig.DEFAULT_NUM_SITES_PER_TYPE;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MarkDupsJitterAnalyserConfig
{
    public final String OutputDir;
    public final String RefGenomeMicrosatelliteFile;
    public final int MinMapQuality;
    public final int MaxSitesPerType;

    public static final String MICROSATELLITE_OUTPUT_DIR = "microsatellite_output_dir";
    private static final String REF_GENOME_MICROSATELLITE_FILE = "ref_genome_microsatellites";
    private static final String MIN_MAP_QUALITY = "microsatellite_min_map_quality";
    private static final String MAX_SITES_PER_TYPE = "max_microsatellite_sites_per_type";

    public MarkDupsJitterAnalyserConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = configBuilder.getValue(MICROSATELLITE_OUTPUT_DIR);
        RefGenomeMicrosatelliteFile = configBuilder.getValue(REF_GENOME_MICROSATELLITE_FILE);
        MinMapQuality = configBuilder.getInteger(MIN_MAP_QUALITY);
        MaxSitesPerType = configBuilder.getInteger(MAX_SITES_PER_TYPE);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(MICROSATELLITE_OUTPUT_DIR, false, "Directory for the output of microsatellite jitter analysis.");
        configBuilder.addPath(REF_GENOME_MICROSATELLITE_FILE, false, "File containing ref genome microsatellite sites");
        configBuilder.addInteger(MIN_MAP_QUALITY, "Minimum mapping quality for an alignment to be used in microsatellite jitter analysis", DEFAULT_MIN_MAPPING_QUALITY);
        configBuilder.addInteger(MAX_SITES_PER_TYPE, "Max number of sites per microsatellite unit / length type", DEFAULT_NUM_SITES_PER_TYPE);
    }
}
