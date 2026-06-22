package com.hartwig.hmftools.orange.report.pdfdata;

import org.jetbrains.annotations.Nullable;

public class CuppaChapterData
{
    public final boolean hasPurpleFail;
    @Nullable
    public final String cuppaSummaryPlotPath;

    public CuppaChapterData(final boolean hasPurpleFail, @Nullable final String cuppaSummaryPlotPath)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.cuppaSummaryPlotPath = cuppaSummaryPlotPath;
    }
}
