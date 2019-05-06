package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.Document;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface ReportChapter {

    @NotNull
    String getName();

    default @Nullable String getPageNumberPrefix() {
        return null;
    }

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

    void render(@NotNull final Document reportDocument);
}
