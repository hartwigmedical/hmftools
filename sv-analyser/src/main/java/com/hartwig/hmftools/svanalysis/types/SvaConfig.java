package com.hartwig.hmftools.svanalysis.types;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SvaConfig
{

    final public int ProximityDistance;
    final public String OutputCsvPath;
    final public String FragileSiteFile;
    final public String LineElementFile;
    final public String ReplicationOriginsFile;
    final public String SampleId;
    final public int MaxSamples;
    final public boolean WriteVisualisationData;
    final public int MaxClusterSize; // for analysis and chaining
    public boolean MergeInconsistentFoldbacks;

    public boolean LogVerbose;
    public String RequiredAnnotations;

    // config options
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String REPLICATION_ORIGINS_FILE = "replication_origins_file";
    private static final String LOG_VERBOSE = "log_verbose";
    private static final String MAX_SAMPLES = "max_samples"; // for testing only
    private static final String WRITE_VISUALISATION_DATA = "write_vis_data";
    private static final String MAX_CLUSTER_SIZE = "max_cluster_size";
    private static final String REQUIRED_ANNOTATIONS = "annotations";
    private static final String LOG_CLUSTER_ID = "log_cluster_id"; // for logging and breakends

    public static int SPECIFIC_CLUSTER_ID = -1;

    private static int DEFAULT_MAX_CLUSTER_SIZE = 1000;

    public SvaConfig(final CommandLine cmd, final String sampleId)
    {
        String dataOutputDir = cmd.getOptionValue(DATA_OUTPUT_PATH);
        if(!dataOutputDir.endsWith(File.separator))
            dataOutputDir += File.separator;

        OutputCsvPath = dataOutputDir;

        ProximityDistance = Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0"));
        SampleId = sampleId;
        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        ReplicationOriginsFile = cmd.getOptionValue(REPLICATION_ORIGINS_FILE, "");
        RequiredAnnotations = cmd.getOptionValue(REQUIRED_ANNOTATIONS, "");
        MaxSamples = Integer.parseInt(cmd.getOptionValue(MAX_SAMPLES, "0"));
        LogVerbose = cmd.hasOption(LOG_VERBOSE);
        WriteVisualisationData = cmd.hasOption(WRITE_VISUALISATION_DATA);
        MaxClusterSize = cmd.hasOption(MAX_CLUSTER_SIZE) ? Integer.parseInt(cmd.getOptionValue(MAX_CLUSTER_SIZE)) : DEFAULT_MAX_CLUSTER_SIZE;
        MergeInconsistentFoldbacks = false;

        SPECIFIC_CLUSTER_ID = Integer.parseInt(cmd.getOptionValue(LOG_CLUSTER_ID, "-1"));
    }

    public SvaConfig(int proximityDistance)
    {
        ProximityDistance = proximityDistance;
        OutputCsvPath = "";
        FragileSiteFile = "";
        LineElementFile = "";
        ReplicationOriginsFile = "";
        RequiredAnnotations = "";
        SampleId = "";
        MaxSamples = 0;
        LogVerbose = false;
        WriteVisualisationData = false;
        MaxClusterSize = DEFAULT_MAX_CLUSTER_SIZE;
        MergeInconsistentFoldbacks = true;
    }

    public boolean hasMultipleSamples() { return SampleId.isEmpty() || SampleId.equals("*"); }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 5000");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file");
        options.addOption(REPLICATION_ORIGINS_FILE, true, "Origins of replication file");
        options.addOption(MAX_SAMPLES, true, "Limit to X samples for testing");
        options.addOption(LOG_VERBOSE, false, "Log extra detail");
        options.addOption(WRITE_VISUALISATION_DATA, false, "Optional: write files for Circos");
        options.addOption(MAX_CLUSTER_SIZE, true, "Optional: max cluster size for chaining");
        options.addOption(REQUIRED_ANNOTATIONS, true, "Optional: string list of annotations");
        options.addOption(LOG_CLUSTER_ID, true, "Optional: log specific cluster details");
    }
}
