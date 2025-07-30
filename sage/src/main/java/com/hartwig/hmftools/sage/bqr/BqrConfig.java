package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MIN_MAP_QUAL;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.SageConstants;

public class BqrConfig
{
    public final boolean Enabled;
    public final boolean LoadBqrFiles;
    public final boolean WriteFile;
    public final boolean WritePlot;
    public final boolean WritePositions;
    public final boolean WriteReads;
    public final int SampleSize;
    public final int MinMapQuality;

    private static final String DISABLE_BQR = "bqr_disable";
    private static final String LOAD_BQR_FILES = "bqr_load";
    private static final String WRITE_BQR_READS = "bqr_write_reads";
    private static final String WRITE_BQR_POSITIONS = "bqr_write_positions";
    private static final String WRITE_BQR_PLOT = "bqr_write_plot";

    private static final String BQR_SAMPLE_SIZE = "bqr_sample_size";
    private static final String BQR_MIN_MAP_QUAL = "bqr_min_map_qual";

    public BqrConfig(final ConfigBuilder configBuilder)
    {
        Enabled = !configBuilder.hasFlag(DISABLE_BQR);

        if(configBuilder.hasFlag(LOAD_BQR_FILES))
        {
            LoadBqrFiles = true;
            WriteFile = false;
            WritePlot = false;
            WritePositions = false;
            WriteReads = false;
        }
        else
        {
            LoadBqrFiles = false;
            WriteFile = Enabled; // written by default
            WritePlot = configBuilder.hasFlag(WRITE_BQR_PLOT);
            WritePositions = configBuilder.hasFlag(WRITE_BQR_POSITIONS);
            WriteReads = configBuilder.hasFlag(WRITE_BQR_READS);
        }

        SampleSize = configBuilder.getInteger(BQR_SAMPLE_SIZE);
        MinMapQuality = configBuilder.getInteger(BQR_MIN_MAP_QUAL);
    }

    public BqrConfig()
    {
        Enabled = false;
        WritePlot = false;
        WriteReads = false;
        WritePositions = false;
        LoadBqrFiles = false;
        WriteFile = false;
        SampleSize = SageConstants.BQR_SAMPLE_SIZE;
        MinMapQuality = DEFAULT_BQR_MIN_MAP_QUAL;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_BQR, "Disable Base Quality Recalibration");
        configBuilder.addFlag(WRITE_BQR_PLOT, "Generate BQR plot");
        configBuilder.addFlag(WRITE_BQR_POSITIONS, "Write positional data as contributes to BQR");
        configBuilder.addFlag(WRITE_BQR_READS, "Write detailed read data as contributes to BQR");
        configBuilder.addFlag(LOAD_BQR_FILES, "Attemps to find and load previously-written BQR files");
        configBuilder.addInteger(BQR_SAMPLE_SIZE, "BQR sampling size per autosome", SageConstants.BQR_SAMPLE_SIZE);
        configBuilder.addInteger(BQR_MIN_MAP_QUAL, "BQR min base quality remap qual", DEFAULT_BQR_MIN_MAP_QUAL);
    }
}
