package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.cfreport.PageEventHandler;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.AreaBreak;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.AreaBreakType;

public abstract class ReportChapter {

    public enum ChapterType {
        SummaryChapter,      // Full height sidepanel, full content
        ClosingChapter,      // Full height sidepanel, just ID and report date
        ContentChapter       // Short sidepanel with just ID and report date
    };

    public void render(PageEventHandler eventHandler, Document report) {

        // Setup event handler
        boolean fullSidebar = getChapterType() == ChapterType.SummaryChapter || getChapterType() == ChapterType.ClosingChapter;
        boolean fullContent = getChapterType() == ChapterType.SummaryChapter;
        eventHandler.setSidebarType(fullSidebar, fullContent);

        // Initialize chapter
        startChapter(getName(), getChapterType() != ChapterType.ClosingChapter.SummaryChapter, report);

        // Render chapter specific content
        renderChapterContent(report);

    }

    private static final void startChapter(String chapterName, boolean addPageBreak, Document report) {

        // Add page break for all subsequent chapters
        if (addPageBreak) {
            report.add(new AreaBreak(AreaBreakType.NEXT_PAGE));
        }

        // Add chapter title
        report.add(new Paragraph(chapterName).addStyle(ReportResources.chapterTitleStyle()));

    }

    public abstract String getName();

    public abstract ChapterType getChapterType();


    protected abstract void renderChapterContent(Document report);


}
