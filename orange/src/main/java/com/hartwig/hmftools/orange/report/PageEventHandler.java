package com.hartwig.hmftools.orange.report;

import com.google.common.io.Resources;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.components.Footer;
import com.hartwig.hmftools.orange.report.components.Header;
import com.hartwig.hmftools.orange.report.components.SidePanel;
import com.itextpdf.kernel.events.Event;
import com.itextpdf.kernel.events.IEventHandler;
import com.itextpdf.kernel.events.PdfDocumentEvent;
import com.itextpdf.kernel.pdf.PdfDocument;
import com.itextpdf.kernel.pdf.PdfOutline;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.navigation.PdfExplicitRemoteGoToDestination;

import org.jetbrains.annotations.NotNull;

public class PageEventHandler implements IEventHandler {

    @NotNull
    private final OrangeReport report;

    @NotNull
    private final Footer footer;
    @NotNull
    private final Header header;

    private String chapterTitle = "Undefined";
    private boolean firstPageOfChapter = true;

    private PdfOutline outline = null;

    PageEventHandler(@NotNull final OrangeReport report) {
        this.report = report;

        this.header = new Header(Resources.getResource("orange_circos.png").getPath());
        this.footer = new Footer();
    }

    @Override
    public void handleEvent(@NotNull Event event) {
        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if (documentEvent.getType().equals(PdfDocumentEvent.START_PAGE)) {
            PdfPage page = documentEvent.getPage();

            header.renderHeader(page);
            if (firstPageOfChapter) {
                firstPageOfChapter = false;

                createChapterBookmark(documentEvent.getDocument(), chapterTitle);
            }

            SidePanel.renderSidePanel(page, report);
            footer.renderFooter(page);
        }
    }

    void chapterTitle(@NotNull String chapterTitle) {
        this.chapterTitle = chapterTitle;
    }

    void resetChapterPageCounter() {
        firstPageOfChapter = true;
    }

    void writeDynamicTextParts(@NotNull PdfDocument document) {
        footer.writeTotalPageCount(document);
    }

    private void createChapterBookmark(@NotNull PdfDocument pdf, @NotNull String title) {
        if (outline == null) {
            outline = pdf.getOutlines(false);
        }

        PdfOutline chapterItem = outline.addOutline(title);
        chapterItem.addDestination(PdfExplicitRemoteGoToDestination.createFitH(pdf.getNumberOfPages(), 0));
    }
}
