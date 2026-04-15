package com.hartwig.hmftools.orange.report.chapters;

import java.io.IOException;

import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public interface ReportChapter
{
    @NotNull
    String name();

    @NotNull
    PDRectangle pageSize();

    default float contentWidth()
    {
        return pageSize().getWidth() - (5 + ReportResources.PAGE_MARGIN_LEFT + ReportResources.PAGE_MARGIN_RIGHT);
    }

    void render(@NotNull DocumentContext document) throws IOException;
}
