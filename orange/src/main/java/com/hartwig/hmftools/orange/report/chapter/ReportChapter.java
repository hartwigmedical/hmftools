package com.hartwig.hmftools.orange.report.chapter;

import com.itextpdf.layout.Document;

import org.jetbrains.annotations.NotNull;

public interface ReportChapter {

    @NotNull
    String name();

    void render(@NotNull Document reportDocument);

}
