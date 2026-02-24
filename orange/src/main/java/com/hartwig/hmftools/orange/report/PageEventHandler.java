package com.hartwig.hmftools.orange.report;

import com.google.common.io.Resources;
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

public class PageEventHandler implements IEventHandler
{
    private final Header mHeader;
    private final Footer mFooter;
    private final SidePanel mSidePanel;

    private String mChapterTitle;
    private boolean mFirstPageOfChapter;
    private PdfOutline mOutline;

    static PageEventHandler create(
            final String sampleId, final String pipelineVersion, final ReportResources reportResources, boolean addDisclaimer)
    {
        return new PageEventHandler(new Header(Resources.getResource("orange_circos.png"), reportResources, addDisclaimer),
                new Footer(reportResources, addDisclaimer),
                new SidePanel(sampleId, pipelineVersion, reportResources));
    }

    private PageEventHandler(final Header header, final Footer footer, final SidePanel sidePanel)
    {
        mHeader = header;
        mFooter = footer;
        mSidePanel = sidePanel;

        mChapterTitle = "Undefined";
        mFirstPageOfChapter = true;
        mOutline = null;
    }

    @Override
    public void handleEvent(final Event event)
    {
        PdfDocumentEvent documentEvent = (PdfDocumentEvent) event;
        if(documentEvent.getType().equals(PdfDocumentEvent.START_PAGE))
        {
            PdfPage page = documentEvent.getPage();

            mHeader.renderHeader(page);
            if(mFirstPageOfChapter)
            {
                mFirstPageOfChapter = false;

                createChapterBookmark(documentEvent.getDocument(), mChapterTitle);
            }

            mSidePanel.renderSidePanel(page);
            mFooter.renderFooter(page);
        }
    }

    void chapterTitle(final String chapterTitle)
    {
        this.mChapterTitle = chapterTitle;
    }

    void resetChapterPageCounter()
    {
        mFirstPageOfChapter = true;
    }

    void writeFooters(final PdfDocument document)
    {
        mFooter.writeFooters(document);
    }

    private void createChapterBookmark(final PdfDocument pdf, final String title)
    {
        if(mOutline == null)
        {
            mOutline = pdf.getOutlines(false);
        }

        PdfOutline chapterItem = mOutline.addOutline(title);
        chapterItem.addDestination(PdfExplicitRemoteGoToDestination.createFitH(pdf.getNumberOfPages(), 0));
    }
}
