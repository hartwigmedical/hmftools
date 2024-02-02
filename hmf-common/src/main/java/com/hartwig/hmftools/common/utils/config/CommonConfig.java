package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public final class CommonConfig
{
    public static final String PLURALS_DESC = "(s), separated by ','";

    public static final String SAMPLE = "sample";
    public static final String SAMPLE_DESC = "Sample ID";

    public static final String TUMOR = "tumor";
    public static final String TUMOR_DESC = "Tumor ID";
    public static final String TUMOR_IDS_DESC = TUMOR_DESC + PLURALS_DESC;

    public static final String REFERENCE = "reference";
    public static final String REFERENCE_DESC = "Reference ID";
    public static final String REFERENCE_IDS_DESC = REFERENCE_DESC + PLURALS_DESC;

    public static final String TUMOR_BAM = "tumor_bam";
    public static final String TUMOR_BAM_DESC = "Tumor BAM file";
    public static final String TUMOR_BAMS_DESC = TUMOR_BAM_DESC + PLURALS_DESC;

    public static final String REFERENCE_BAM = "reference_bam";
    public static final String REFERENCE_BAM_DESC = "Reference BAM file";
    public static final String REFERENCE_BAMS_DESC = REFERENCE_BAM_DESC + PLURALS_DESC;

    public static final String RNA_SAMPLE_ID = "rna_sample";
    public static final String RNA_SAMPLE_DESC = "RNA sample ID";

    public static final String RNA_BAM = "rna_bam";
    public static final String RNA_BAM_DESC = "RNA BAM file";


    public static final String LOG_READ_IDS = "log_read_ids";
    public static final String LOG_READ_IDS_DESC = "Log specific read IDs, separated by ';'";

    public static final String PERF_DEBUG = "perf_debug";
    public static final String PERF_DEBUG_DESC = "Detailed performance tracking and logging";

    public static final String TARGET_REGIONS_BED = "target_regions_bed";
    public static final String TARGET_REGIONS_BED_DESC = "Target regions BED file";

    public static final String SAMPLE_DATA_DIR_CFG = "sample_data_dir";
    public static final String SAMPLE_DATA_DIR_DESC = "Path to sample pipeline files";

    public static final String PIPELINE_SAMPLE_ROOT_DIR = "pipeline_sample_root_dir";
    public static final String PIPELINE_SAMPLE_ROOT_DESC = "Path to pipeline sample root directory, expecting sub-directories per tool";

    public static final String AMBER_DIR_CFG = toolDirectory("amber");
    public static final String AMBER_DIR_DESC = toolDirectoryDesc("Amber");

    public static final String CHORD_DIR_CFG = toolDirectory("chord");
    public static final String CHORD_DIR_DESC = toolDirectoryDesc("Chord");

    public static final String COBALT_DIR_CFG = toolDirectory("cobalt");
    public static final String COBALT_DIR_DESC = toolDirectoryDesc("Cobalt");

    public static final String CUPPA_DIR_CFG = toolDirectory("cuppa");
    public static final String CUPPA_DIR_DESC = toolDirectoryDesc("Cuppa");

    public static final String ISOFOX_DIR_CFG = toolDirectory("isofox");
    public static final String ISOFOX_DIR_DESC = toolDirectoryDesc("Isofox");

    public static final String LILAC_DIR_CFG = toolDirectory("lilac");
    public static final String LILAC_DIR_DESC = toolDirectoryDesc("Lilac");

    public static final String LINX_DIR_CFG = toolDirectory("linx");
    public static final String LINX_DIR_DESC = toolDirectoryDesc("Linx");

    public static final String LINX_PLOT_DIR_CFG = toolDirectory("linx_plot");
    public static final String LINX_PLOT_DIR_DESC = toolPlotsDirectoryDesc("Linx");

    public static final String LINX_GERMLINE_DIR_CFG = toolDirectory("linx_germline");
    public static final String LINX_GERMLINE_DIR_DESC = toolDirectoryDesc("Linx germline");

    public static final String NEO_DIR_CFG = toolDirectory("neo");
    public static final String NEO_DIR_DESC = toolDirectoryDesc("Neo");

    public static final String PEACH_DIR_CFG = toolDirectory("peach");
    public static final String PEACH_DIR_DESC = toolDirectoryDesc("Peach");

    public static final String PURPLE_DIR_CFG = toolDirectory("purple");
    public static final String PURPLE_DIR_DESC = toolDirectoryDesc("Purple");

    public static final String PURPLE_PLOT_DIR_CFG = toolDirectory("purple_plot");
    public static final String PURPLE_PLOT_DIR_DESC = toolPlotsDirectoryDesc("Purple");

    public static final String SAGE_DIR_CFG = toolDirectory("sage");
    public static final String SAGE_DIR_DESC = toolDirectoryDesc("Sage");

    public static final String SAGE_GERMLINE_DIR_CFG = toolDirectory("sage_germline");
    public static final String SAGE_GERMLINE_DIR_DESC = toolDirectoryDesc("Sage germline");

    public static final String SIGS_DIR_CFG = toolDirectory("sigs");
    public static final String SIGS_DIR_DESC = toolDirectoryDesc("Signatures");

    public static final String VIRUS_DIR_CFG = toolDirectory("virus");
    public static final String VIRUS_DIR_DESC = toolDirectoryDesc("Virus");

    private static String toolDirectory(final String toolName) { return format("%s_dir", toolName); }
    private static String toolDirectoryDesc(final String toolName) { return format("Path to %s pipeline files", toolName); }
    private static String toolPlotsDirectoryDesc(final String toolName) { return format("Path to %s plots", toolName); }

    public static List<String> parseLogReadIds(final ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(LOG_READ_IDS))
            return Arrays.stream(configBuilder.getValue(LOG_READ_IDS).split(ITEM_DELIM, -1)).collect(Collectors.toList());
        else
            return Collections.emptyList();

    }
}
