package com.hartwig.hmftools.svanalysis.analysis;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SvaConfig
{

    final public int ProximityDistance;
    final public String OutputCsvPath;
    final public String FragileSiteFile;
    final public String LineElementFile;
    final public String ReplicationOriginsFile;
    final public String LOHDataFile;
    final public String SampleId;
    final public int MaxSamples;
    final public boolean WriteVisualisationData;

    public boolean LogVerbose;

    // config options
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String REPLICATION_ORIGINS_FILE = "replication_origins_file";
    private static final String LOH_DATA_FILE = "loh_file";
    private static final String LOG_VERBOSE = "log_verbose";
    private static final String MAX_SAMPLES = "max_samples"; // for testing only
    private static final String WRITE_VISUALISATION_DATA = "write_vis_data";

    public SvaConfig(final CommandLine cmd, final String sampleId)
    {
        OutputCsvPath = cmd.getOptionValue(DATA_OUTPUT_PATH);
        ProximityDistance = Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0"));
        SampleId = sampleId;
        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        ReplicationOriginsFile = cmd.getOptionValue(REPLICATION_ORIGINS_FILE, "");
        LOHDataFile = cmd.getOptionValue(LOH_DATA_FILE, "");
        MaxSamples = Integer.parseInt(cmd.getOptionValue(MAX_SAMPLES, "0"));
        LogVerbose = cmd.hasOption(LOG_VERBOSE);
        WriteVisualisationData = cmd.hasOption(WRITE_VISUALISATION_DATA);
    }

    public SvaConfig(int proximityDistance)
    {
        ProximityDistance = proximityDistance;
        OutputCsvPath = "";
        FragileSiteFile = "";
        LineElementFile = "";
        ReplicationOriginsFile = "";
        LOHDataFile = "";
        SampleId = "";
        MaxSamples = 0;
        LogVerbose = false;
        WriteVisualisationData = false;
    }

    public boolean hasMultipleSamples() { return SampleId.isEmpty() || SampleId.equals("*"); }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 5000");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file");
        options.addOption(REPLICATION_ORIGINS_FILE, true, "Origins of replication file");
        options.addOption(LOH_DATA_FILE, true, "Copy Number LOH data file");
        options.addOption(MAX_SAMPLES, true, "Limit to X samples for testing");
        options.addOption(LOG_VERBOSE, false, "Log extra detail");
        options.addOption(WRITE_VISUALISATION_DATA, false, "Optional: write files for Circos");
    }
}
