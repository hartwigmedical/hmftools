package com.hartwig.hmftools.svanalysis.analysis;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SvaConfig
{

    final public int ProximityDistance;
    final public String OutputCsvPath;
    final public String SvPONFile;
    final public String FragileSiteFile;
    final public String LineElementFile;
    final public String LOHDataFile;
    final public String SampleId;

    public boolean LogVerbose;

    // config options
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String SV_PON_FILE = "sv_pon_file";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String LOH_DATA_FILE = "loh_file";
    private static final String LOG_VERBOSE = "log_verbose";

    public SvaConfig(final CommandLine cmd, final String sampleId)
    {
        OutputCsvPath = cmd.getOptionValue(DATA_OUTPUT_PATH);
        ProximityDistance = Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0"));
        SampleId = sampleId;
        SvPONFile = cmd.getOptionValue(SV_PON_FILE, "");
        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        LOHDataFile = cmd.getOptionValue(LOH_DATA_FILE, "");
        LogVerbose = cmd.hasOption(LOG_VERBOSE);
    }

    public SvaConfig(int proximityDistance)
    {
        ProximityDistance = proximityDistance;
        OutputCsvPath = "";
        SvPONFile = "";
        FragileSiteFile = "";
        LineElementFile = "";
        LOHDataFile = "";
        SampleId = "";
        LogVerbose = false;
    }

    public boolean hasMultipleSamples() { return SampleId.isEmpty() || SampleId.equals("*"); }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 5000");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file for SVs");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file for SVs");
        options.addOption(LOH_DATA_FILE, true, "Copy Number LOH data file");
        options.addOption(LOG_VERBOSE, false, "Log extra detail");
    }
}
