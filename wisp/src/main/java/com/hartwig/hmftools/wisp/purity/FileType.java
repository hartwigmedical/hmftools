package com.hartwig.hmftools.wisp.purity;

public enum FileType
{
    SUMMARY("summary", false),
    SOMATICS("somatic_variants", false),
    CN_SEGMENT("cn_segments", false),
    SOMATIC_PEAK("somatic_peak", true),
    CN_PLOT_CALCS("cn_plot_calcs", true);

    private final String mFileId;
    private final boolean mIsPlotData;

    FileType(final String fileId, boolean isPlotData)
    {
        mFileId = fileId;
        mIsPlotData = isPlotData;
    }

    public String fileId() { return mFileId; }
    public boolean isPlotData() { return mIsPlotData; }
}
