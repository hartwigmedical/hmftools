package com.hartwig.hmftools.orange.report.pdfdata;

import org.jetbrains.annotations.Nullable;

public class QualityControlData
{
    public final boolean hasPurpleFail;
    @Nullable
    public final String qSeePlotPath;

    public QualityControlData(final boolean hasPurpleFail, @Nullable final String qSeePlotPath)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.qSeePlotPath = qSeePlotPath;
    }
}
