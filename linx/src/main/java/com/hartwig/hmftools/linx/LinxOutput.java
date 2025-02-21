package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.linx.WriteType.CLUSTERING;
import static com.hartwig.hmftools.linx.WriteType.LINK;
import static com.hartwig.hmftools.linx.WriteType.STANDARD_TYPES;
import static com.hartwig.hmftools.linx.WriteType.SV_DATA;
import static com.hartwig.hmftools.linx.WriteType.VIS_DATA;
import static com.hartwig.hmftools.linx.WriteType.WRITE_STANDARD;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class LinxOutput
{
    public final List<WriteType> WriteTypes;
    public final boolean WriteCohortFiles;
    public final boolean WriteSingleSVClusters;

    public final int LogChainingMaxSize;

    private final boolean mWriteAll;

    private static final String WRITE_TYPES = "write_type";

    private static final String WRITE_ALL = "write_all";
    private static final String WRITE_SV_DATA = "write_sv_data";
    private static final String WRITE_SINGLE_SV_CLUSTERS = "write_single_sv_clusters";
    private static final String WRITE_CLUSTER_HISTORY = "write_cluster_history";
    private static final String WRITE_LINKS = "write_links";
    private static final String WRITE_VIS_DATA = "write_vis_data";
    private static final String WRITE_COHORT_FILES = "write_cohort";
    private static final String NO_VIS_FILES = "no_vis_files";
    private static final String LOG_CHAIN_MAX_SIZE = "log_chain_size";

    public static final char ITEM_DELIM_CHR = ITEM_DELIM.charAt(0);

    public LinxOutput(final ConfigBuilder configBuilder, boolean defaultWrite)
    {
        WriteTypes = WriteType.parseConfig(configBuilder.getValue(WRITE_TYPES));

        mWriteAll = configBuilder.getValue(WRITE_TYPES).equals(WRITE_ALL);

        if(mWriteAll)
        {
            WriteSingleSVClusters = true;
        }
        else
        {
            WriteSingleSVClusters = configBuilder.hasFlag(WRITE_SINGLE_SV_CLUSTERS) || defaultWrite;
        }

        WriteCohortFiles = configBuilder.hasFlag(WRITE_COHORT_FILES);

        LogChainingMaxSize = configBuilder.getInteger(LOG_CHAIN_MAX_SIZE);
    }

    public boolean writeAll() { return mWriteAll; }
    public boolean writeVisualisationData() { return WriteTypes.contains(VIS_DATA); }
    public boolean writeSvData() { return WriteTypes.contains(SV_DATA); }
    public boolean writeLinks() { return WriteTypes.contains(LINK); }
    public boolean writeClusterHistory() { return WriteTypes.contains(CLUSTERING); }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(
                WRITE_TYPES, false, enumValueSelectionAsStr(WriteType.values(), "Write types"), WRITE_STANDARD);

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
        mWriteAll = false;
        WriteTypes = STANDARD_TYPES;
        WriteSingleSVClusters = false;
        WriteCohortFiles = false;
        LogChainingMaxSize = 0;
    }
}
