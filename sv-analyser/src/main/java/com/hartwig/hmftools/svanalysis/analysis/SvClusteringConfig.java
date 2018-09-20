package com.hartwig.hmftools.svanalysis.analysis;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SvClusteringConfig {

    final public int ClusterBaseDistance;
    final public String OutputCsvPath;
    final public String SvPONFile;
    final public String FragileSiteFile;
    final public String LineElementFile;
    final public String ExternalAnnotationsFile;
    final public String GeneDataFile;
    final public boolean UseCombinedOutputFile;
    final public boolean UseGridss;

    // config options
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String SV_PON_FILE = "sv_pon_file";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";
    private static final String EXTERNAL_SV_DATA_FILE = "ext_sv_data_file";
    private static final String EXTERNAL_DATA_LINK_FILE = "ext_data_link_file";
    private static final String DRIVER_GENES_FILE = "driver_gene_file";
    private static final String USE_GRIDSS = "use_gridss";


    public static final int DEFAULT_BASE_DISTANCE = 100000;

    public SvClusteringConfig(final CommandLine cmd, final String sampleId)
    {
        OutputCsvPath = cmd.getOptionValue(DATA_OUTPUT_PATH);
        ClusterBaseDistance = Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0"));
        UseCombinedOutputFile = sampleId.isEmpty() || sampleId.equals("*");
        SvPONFile = cmd.getOptionValue(SV_PON_FILE, "");
        FragileSiteFile = cmd.getOptionValue(FRAGILE_SITE_FILE, "");
        LineElementFile = cmd.getOptionValue(LINE_ELEMENT_FILE, "");
        ExternalAnnotationsFile = cmd.getOptionValue(EXTERNAL_SV_DATA_FILE, "");
        GeneDataFile = cmd.getOptionValue(DRIVER_GENES_FILE, "");
        UseGridss = cmd.hasOption(USE_GRIDSS);
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 1000");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file for SVs");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file for SVs");
        options.addOption(EXTERNAL_SV_DATA_FILE, true, "External file with per-SV annotations");
        options.addOption(EXTERNAL_DATA_LINK_FILE, true, "External SV data file, mapped by position info");
        options.addOption(DRIVER_GENES_FILE, true, "Gene data file");
        options.addOption(USE_GRIDSS, false, "Using GRIDSS variant calling input data");
    }
}
