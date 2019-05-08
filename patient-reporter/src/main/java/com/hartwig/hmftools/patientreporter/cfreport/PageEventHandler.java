package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.PatientReport;
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

class PageEventHandler implements IEventHandler {

    @NotNull
    private final SampleReport sampleReport;
    @NotNull
    private final Footer footer;
    @NotNull
    private final Header header;

    private boolean fullSidebar;
    private boolean fullSidebarContent;

    private String chapterTitle = "Undefined";
    private String pageNumberPrefix = null;
    private boolean firstPageOfChapter = true;

    private PdfOutline outline = null;

    PageEventHandler(@NotNull final PatientReport patientReport) {
        this.sampleReport = patientReport.sampleReport();
        this.header = new Header(patientReport.logoCompanyPath());
        this.footer = new Footer();
    }

    @Override
    public void handleEvent(@NotNull Event event) {
        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType().equals(PdfDocumentEvent.START_PAGE)) {
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

    public void chapterTitle(@NotNull String chapterTitle) {
        this.chapterTitle = chapterTitle;
    }

    public void pageNumberPrefix(@Nullable String pageNumberPrefix) {
        this.pageNumberPrefix = pageNumberPrefix;
    }

    public void sidebarType(boolean full, boolean fullContent) {
        fullSidebar = full;
        fullSidebarContent = fullSidebar && fullContent;
    }

    public void resetChapterPageCounter() {
        firstPageOfChapter = true;
    }

    public void writeDynamicTextParts(@NotNull PdfDocument document) {
        header.writeChapterTitles(document);
        footer.writeTotalPageCount(document);
    }

    private void createChapterBookmark(@NotNull PdfDocument pdf, @NotNull String title) {
        if (outline == null) {
            outline = pdf.getOutlines(false);
        }

        final PdfOutline chapterItem = outline.addOutline(title);
        chapterItem.addDestination(PdfExplicitRemoteGoToDestination.createFitH(pdf.getNumberOfPages(), 0));
    }
}
