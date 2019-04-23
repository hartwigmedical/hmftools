package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.Document;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public interface ReportChapter {

    String getName();

    default boolean isFullWidth() {
        return true;
    }

    default boolean hasCompleteSidebar() {
        return false;
    }

    default float getContentWidth() {
        return isFullWidth()
                ? ReportResources.CONTENT_WIDTH_WIDE
                : ReportResources.CONTENT_WIDTH_NARROW;
    }

    void render(@NotNull final Document reportDocument) throws IOException;

}
