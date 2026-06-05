package com.hartwig.hmftools.orange.report;

import com.google.common.io.Resources;
import com.hartwig.hmftools.orange.report.components.Footer;
import com.hartwig.hmftools.orange.report.components.Header;
import com.hartwig.hmftools.orange.report.components.SidePanel;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.interactive.documentnavigation.destination.PDPageFitWidthDestination;
import org.apache.pdfbox.pdmodel.interactive.documentnavigation.outline.PDDocumentOutline;
import org.apache.pdfbox.pdmodel.interactive.documentnavigation.outline.PDOutlineItem;

public class PageEventHandler
{
    private final Header mHeader;
    private final Footer mFooter;
    private final SidePanel mSidePanel;

    private String mChapterTitle;
    private boolean mFirstPageOfChapter;
    private PDDocumentOutline mOutline;

    static PageEventHandler create(final String sampleId, final ReportResources reportResources, boolean addDisclaimer,
            final PDDocument document)
    {
        return new PageEventHandler(
                new Header(Resources.getResource("orange_circos.png"), reportResources, addDisclaimer, document),
                new Footer(reportResources, addDisclaimer),
                new SidePanel(sampleId, reportResources));
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

    void onPageStart(final PDPage page, final PDDocument document)
    {
        mHeader.renderHeader(page, document);

        if(mFirstPageOfChapter)
        {
            mFirstPageOfChapter = false;
            createChapterBookmark(document, page, mChapterTitle);
        }

        mSidePanel.renderSidePanel(page, document);
        mFooter.renderFooter(page, document);
    }

    void chapterTitle(final String chapterTitle)
    {
        this.mChapterTitle = chapterTitle;
    }

    void resetChapterPageCounter()
    {
        mFirstPageOfChapter = true;
    }

    void writeFooters(final PDDocument document)
    {
        mFooter.writeFooters(document);
    }

    private void createChapterBookmark(final PDDocument document, final PDPage page, final String title)
    {
        if(mOutline == null)
        {
            mOutline = new PDDocumentOutline();
            document.getDocumentCatalog().setDocumentOutline(mOutline);
        }

        PDPageFitWidthDestination dest = new PDPageFitWidthDestination();
        dest.setPage(page);

        PDOutlineItem chapterItem = new PDOutlineItem();
        chapterItem.setTitle(title);
        chapterItem.setDestination(dest);
        mOutline.addLast(chapterItem);
    }
}
