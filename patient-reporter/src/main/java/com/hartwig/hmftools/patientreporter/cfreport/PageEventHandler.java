package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfPage;
import org.jetbrains.annotations.NotNull;

public class PageEventHandler implements IEventHandler {

    private SampleReport sampleReport;

    private Footer footer;
    private Header header;

    private boolean fullSidebar;
    private boolean fullSidebarContent;

    private String chapterTitle = "Undefined";
    private boolean firstPageOfChapter = true;

    public PageEventHandler(@NotNull final SampleReport sampleReport) {
        this.sampleReport = sampleReport;
        this.header = new Header();
        this.footer = new Footer();
    }

    public void setChapterTitle(String chapterTitle) {
        this.chapterTitle = chapterTitle;
    }

    public void setSidebarType(boolean full, boolean fullContent) {
        fullSidebar = full;
        fullSidebarContent = fullSidebar ? fullContent : false;
    }

    public void resetChapterPageCounter() {
        firstPageOfChapter = true;
    }

    public void writeTotalPageCount(@NotNull PdfDocument document) {
        footer.writeTotalPageCount(document);
    }

    @Override
    public void handleEvent(Event event) {

        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType() == PdfDocumentEvent.START_PAGE) {

            final PdfPage page = documentEvent.getPage();

            header.renderHeader(chapterTitle, firstPageOfChapter, page);
            if (firstPageOfChapter) {
                firstPageOfChapter = false;
            }
            SidePanel.renderSidePanel(page, sampleReport, fullSidebar, fullSidebarContent);
            footer.renderFooter(page, !fullSidebar);

        }

    }

}
