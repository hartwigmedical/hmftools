package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_BQR_MIN_MAP_QUAL;

import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.SageConstants;

public class BqrConfig
{
    public final boolean Enabled;
    public final boolean LoadBqrFiles;
    public final boolean UsePanel;
    public final boolean WriteFile;
    public final boolean WritePlot;
    public final boolean WritePositions;
    public final boolean WriteReads;
    public final boolean FullBam;
    public final boolean ExcludeKnown;
    public final int SampleSize;
    public final int MinMapQuality;

    private static final String DISABLE_BQR = "bqr_disable";
    private static final String LOAD_BQR_FILES = "bqr_load";
    private static final String FULL_BAM = "bqr_full_bam";
    private static final String EXCLUDE_KNOWN_VARIANTS = "bqr_exclude_known";
    private static final String USE_PANEL = "bqr_use_panel";
    private static final String WRITE_BQR_READS = "bqr_write_reads";
    private static final String WRITE_BQR_POSITIONS = "bqr_write_positions";
    private static final String WRITE_BQR_PLOT = "bqr_write_plot";

    @Deprecated
    private static final String WRITE_BQR_PLOT_OLD = "write_bqr_plot";

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
            WritePlot = configBuilder.hasFlag(WRITE_BQR_PLOT) || configBuilder.hasFlag(WRITE_BQR_PLOT_OLD);
            WritePositions = configBuilder.hasFlag(WRITE_BQR_POSITIONS);
            WriteReads = configBuilder.hasFlag(WRITE_BQR_READS);
        }

        FullBam = configBuilder.hasFlag(FULL_BAM);
        ExcludeKnown = configBuilder.hasFlag(EXCLUDE_KNOWN_VARIANTS);
        UsePanel = configBuilder.hasFlag(USE_PANEL);

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
        ExcludeKnown = false;
        FullBam = false;
        UsePanel = false;
    }

    public static boolean useReadType(final SageConfig config)
    {
        return config.Quality.HighDepthMode || config.Sequencing.HasUMIs || config.Sequencing.Type == SequencingType.ULTIMA;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(DISABLE_BQR, "Disable Base Quality Recalibration");
        configBuilder.addFlag(WRITE_BQR_PLOT, "Generate BQR plot");
        configBuilder.addFlag(WRITE_BQR_PLOT_OLD, "Decpreated, use instead: " + WRITE_BQR_PLOT);
        configBuilder.addFlag(WRITE_BQR_POSITIONS, "Write positional data as contributes to BQR");
        configBuilder.addFlag(WRITE_BQR_READS, "Write detailed read data as contributes to BQR");
        configBuilder.addFlag(LOAD_BQR_FILES, "Attemps to find and load previously-written BQR files");
        configBuilder.addFlag(FULL_BAM, "Run over full BAM");
        configBuilder.addFlag(USE_PANEL, "Run on panel only");
        configBuilder.addFlag(EXCLUDE_KNOWN_VARIANTS, "Exclude known variants in append mode");
        configBuilder.addInteger(BQR_SAMPLE_SIZE, "BQR sampling size per autosome", SageConstants.BQR_SAMPLE_SIZE);
        configBuilder.addInteger(BQR_MIN_MAP_QUAL, "BQR min base quality remap qual", DEFAULT_BQR_MIN_MAP_QUAL);
    }
}
