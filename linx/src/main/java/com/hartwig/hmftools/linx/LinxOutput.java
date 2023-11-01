package com.hartwig.hmftools.linx;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class LinxOutput
{
    public final boolean WriteAll;
    public final boolean WriteSvData; // all SV table fields to cohort file
    public final boolean WriteVisualisationData;

    public final boolean WriteClusterHistory;
    public final boolean WriteSingleSVClusters;
    public final boolean WriteLinks;
    public final boolean WriteCohortFiles;

    public final int LogChainingMaxSize;

    private static final String WRITE_ALL = "write_all";
    private static final String WRITE_SV_DATA = "write_sv_data";
    private static final String WRITE_SINGLE_SV_CLUSTERS = "write_single_sv_clusters";
    private static final String WRITE_CLUSTER_HISTORY = "write_cluster_history";
    private static final String WRITE_LINKS = "write_links";
    private static final String WRITE_VIS_DATA = "write_vis_data";
    private static final String WRITE_COHORT_FILES = "write_cohort";
    private static final String NO_VIS_FILES = "no_vis_files";
    private static final String LOG_CHAIN_MAX_SIZE = "log_chain_size";

    public static final char ITEM_DELIM_CHR = ';';

    public LinxOutput(final ConfigBuilder configBuilder, boolean defaultWrite)
    {
        WriteAll = configBuilder.hasFlag(WRITE_ALL);

        if(WriteAll)
        {
            WriteSvData = true;
            WriteVisualisationData = true;
            WriteClusterHistory = true;
            WriteSingleSVClusters = true;
            WriteLinks = true;
        }
        else
        {
            WriteVisualisationData = (configBuilder.hasFlag(WRITE_VIS_DATA) || defaultWrite) && !configBuilder.hasFlag(NO_VIS_FILES);

            WriteSvData = configBuilder.hasFlag(WRITE_SV_DATA) || defaultWrite;
            WriteClusterHistory = configBuilder.hasFlag(WRITE_CLUSTER_HISTORY);
            WriteSingleSVClusters = configBuilder.hasFlag(WRITE_SINGLE_SV_CLUSTERS) || defaultWrite;
            WriteLinks = configBuilder.hasFlag(WRITE_LINKS) || defaultWrite;
        }

        WriteCohortFiles = configBuilder.hasFlag(WRITE_COHORT_FILES);

        LogChainingMaxSize = configBuilder.getInteger(LOG_CHAIN_MAX_SIZE);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addFlag(WRITE_ALL, "Optional: write all batch-run output files");
        configBuilder.addFlag(WRITE_SV_DATA, "Optional: include all SV table fields (batch-mode)");
        configBuilder.addFlag(WRITE_LINKS, "Optional: write chain links (batch-mode)");
        configBuilder.addFlag(WRITE_CLUSTER_HISTORY, "Optional: write clustering history (batch-mode)");
        configBuilder.addFlag(WRITE_SINGLE_SV_CLUSTERS, "Optional: write cluster data for single SV clusters (batch-mode)");
        configBuilder.addFlag(WRITE_VIS_DATA, "Optional: write files for Circos (batch-mode)");
        configBuilder.addFlag(WRITE_COHORT_FILES, "Optional: write cohort files even for single sample");
        configBuilder.addFlag(NO_VIS_FILES, "Disable visualiser files output");

        configBuilder.addInteger(
                LOG_CHAIN_MAX_SIZE,
                "Write file with chaining diagnostics for chains less than this, zero is disabled", 0);
    }

    public LinxOutput()
    {
        WriteAll = false;
        WriteSvData = false;
        WriteVisualisationData = false;
        WriteClusterHistory = false;
        WriteSingleSVClusters = false;
        WriteLinks = false;
        WriteCohortFiles = false;

        LogChainingMaxSize = 0;
    }
}
