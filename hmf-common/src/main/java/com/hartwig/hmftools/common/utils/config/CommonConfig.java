package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class CommonConfig
{
    public static final String SAMPLE = "sample";
    public static final String SAMPLE_DESC = "Sample ID";

    public static final String TUMOR = "tumor";
    public static final String TUMOR_DESC = "Tumor ID";

    public static final String REFERENCE = "reference";
    public static final String REFERENCE_DESC = "Reference ID";

    public static final String LOG_READ_IDS = "log_read_ids";
    public static final String LOG_READ_IDS_DESC = "Log specific read IDs, separated by ';'";

    public static final String PERF_DEBUG = "perf_debug";
    public static final String PERF_DEBUG_DESC = "Detailed performance tracking and logging";

    public static final String PURPLE_DIR_CFG = toolDirectory("purple");
    public static final String PURPLE_DIR_DESC = toolDirectoryDesc("purple");

    public static final String AMBER_DIR_CFG = toolDirectory("amber");
    public static final String AMBER_DIR_DESC = toolDirectoryDesc("amber");

    public static final String COBALT_DIR_CFG = toolDirectory("cobalt");
    public static final String COBALT_DIR_DESC = toolDirectoryDesc("cobalt");

    public static final String LINX_DIR_CFG = toolDirectory("linx");
    public static final String LINX_DIR_DESC = toolDirectoryDesc("linx");

    public static final String CHORD_DIR_CFG = toolDirectory("chord");
    public static final String CHORD_DIR_DESC = toolDirectoryDesc("chord");

    public static final String LILAC_DIR_CFG = toolDirectory("lilac");
    public static final String LILAC_DIR_DESC = toolDirectoryDesc("lilac");

    public static final String ISOFOX_DIR_CFG = toolDirectory("isofox");
    public static final String ISOFOX_DIR_DESC = toolDirectoryDesc("isofox");

    public static final String NEO_DIR_CFG = toolDirectory("neo");
    public static final String NEO_DIR_DESC = toolDirectoryDesc("neo");

    private static String toolDirectory(final String toolName) { return format("%s_dir", toolName); }
    private static String toolDirectoryDesc(final String toolName) { return format("Path to %s pipeline files", toolName); }

    public static List<String> parseLogReadIds(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(LOG_READ_IDS))
            return Arrays.stream(configBuilder.getValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList());
        else
            return Collections.emptyList();

    }
}
