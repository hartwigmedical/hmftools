package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfPage;

public class PageEventHandler implements IEventHandler {

    private boolean fullSidebar;
    private boolean fullSidebarContent;

    private String chapterTitle = "Undefined";
    private boolean firstPageOfChapter = true;

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

    @Override
    public void handleEvent(Event event) {

        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType() == PdfDocumentEvent.START_PAGE) {

            final PdfPage page = documentEvent.getPage();

            Header.addHeader(chapterTitle, firstPageOfChapter, page);
            if (firstPageOfChapter) {
                firstPageOfChapter = false;
            }
            SidePanel.renderSidePanel(page, fullSidebar, fullSidebarContent);
            Footer.renderFooter(page, !fullSidebar);
        }

    }

}
