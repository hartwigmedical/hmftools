package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;

import org.jetbrains.annotations.NotNull;

public interface ReportChapter {

    @NotNull
    String name();

    @NotNull
    PageSize pageSize();

    default float contentWidth() {
        return pageSize().getWidth() - (5 + ReportResources.PAGE_MARGIN_LEFT + ReportResources.PAGE_MARGIN_RIGHT);
    }

    void render(@NotNull Document document);

}
