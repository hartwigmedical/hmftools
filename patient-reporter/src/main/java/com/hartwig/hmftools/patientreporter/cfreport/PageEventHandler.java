package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.cfreport.components.Footer;
import com.hartwig.hmftools.patientreporter.cfreport.components.Header;
import com.hartwig.hmftools.patientreporter.cfreport.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfOutline;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.navigation.PdfExplicitRemoteGoToDestination;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PageEventHandler implements IEventHandler {

    private SampleReport sampleReport;

    private Footer footer;
    private Header header;

    private boolean fullSidebar;
    private boolean fullSidebarContent;

    private String chapterTitle = "Undefined";
    private String pageNumberPrefix = null;
    private boolean firstPageOfChapter = true;

    private PdfOutline outline = null;

    public PageEventHandler(@NotNull final SampleReport sampleReport) {
        this.sampleReport = sampleReport;
        this.header = new Header();
        this.footer = new Footer();
    }

    public void setChapterTitle(@NotNull String chapterTitle) {
        this.chapterTitle = chapterTitle;
    }

    public void setPageNumberPrefix(@Nullable String pageNumberPrefix) {
        this.pageNumberPrefix = pageNumberPrefix;
    }

    public void setSidebarType(boolean full, boolean fullContent) {
        fullSidebar = full;
        fullSidebarContent = fullSidebar ? fullContent : false;
    }

    public void resetChapterPageCounter() {
        firstPageOfChapter = true;
    }

    public void writeDynamicTextParts(@NotNull PdfDocument document) {
        header.writeChapterTitles(document);
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

                createChapterBookmark(documentEvent.getDocument(), chapterTitle);

            }
            SidePanel.renderSidePanel(page, sampleReport, fullSidebar, fullSidebarContent);
            footer.renderFooter(page, !fullSidebar, pageNumberPrefix);
        }
    }

    private void createChapterBookmark(PdfDocument pdf, String title) {
        if (outline == null) {
            outline = pdf.getOutlines(false);
        }

        PdfOutline chapterItem = outline.addOutline(title);
        chapterItem.addDestination(PdfExplicitRemoteGoToDestination.createFitH(pdf.getNumberOfPages(), 0));
    }
}
