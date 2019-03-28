package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.cfreport.PageEventHandler;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.property.AreaBreakType;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public abstract class ReportChapter {

    public enum ChapterType {
        SummaryChapter,      // Full height sidepanel, full content
        ClosingChapter,      // Full height sidepanel, just ID and report date
        ContentChapter       // Short sidepanel with just ID and report date
    }

    public void render(@NotNull final PageEventHandler eventHandler, @NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) throws IOException {

        // Reconfigure event handler for upcoming pages
        eventHandler.setChapterTitle(getName());
        eventHandler.resetChapterPageCounter();

        boolean fullSidebar = getChapterType() == ChapterType.SummaryChapter || getChapterType() == ChapterType.ClosingChapter;
        boolean fullContent = getChapterType() == ChapterType.SummaryChapter;
        eventHandler.setSidebarType(fullSidebar, fullContent);

        // Start chapter
        boolean addPageBreak = (getChapterType() != ChapterType.SummaryChapter);
        startChapter(addPageBreak, reportDocument);

        // Render chapter specific content
        renderChapterContent(patientReport, reportDocument);

    }

    private static void startChapter(boolean addPageBreak, @NotNull Document report) {

        // Add page break for all subsequent chapters
        if (addPageBreak) {
            report.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
        }

    }

    public abstract String getName();

    public abstract ChapterType getChapterType();

    protected abstract void renderChapterContent(@NotNull final AnalysedPatientReport patientReport, @NotNull final Document reportDocument) throws IOException;

    final float getContentWidth() {
        if (getChapterType() == ChapterType.SummaryChapter || getChapterType() == ChapterType.ClosingChapter) {
            return ReportResources.CONTENT_WIDTH_NARROW;
        } else {
            return ReportResources.CONTENT_WIDTH_WIDE;
        }
    }

}
