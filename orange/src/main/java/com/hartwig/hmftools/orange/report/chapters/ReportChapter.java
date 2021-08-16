package com.hartwig.hmftools.orange.report.chapters;

import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;

import org.jetbrains.annotations.NotNull;

public interface ReportChapter {

    @NotNull
    String name();

    @NotNull
    PageSize pageSize();

    void render(@NotNull Document document);

}
