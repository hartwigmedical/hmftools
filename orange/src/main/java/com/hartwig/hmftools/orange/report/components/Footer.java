package com.hartwig.hmftools.orange.report.components;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;

public class Footer
{
    private final List<Integer> mPageNumbers;
    private final ReportResources mReportResources;
    private final boolean mAddDisclaimer;

    public Footer(final ReportResources reportResources, boolean addDisclaimer)
    {
        mPageNumbers = Lists.newArrayList();
        mReportResources = reportResources;
        mAddDisclaimer = addDisclaimer;
    }

    /**
     * Record that a page was created (to be numbered later).
     */
    public void renderFooter(final PDPage page, final PDDocument document)
    {
        int pageNumber = document.getPages().indexOf(page) + 1;
        mPageNumbers.add(pageNumber);
    }

    /**
     * Write page numbers on all pages (second pass, once total count is known).
     */
    public void writeFooters(final PDDocument document)
    {
        int totalPageCount = document.getNumberOfPages();
        for(int i = 0; i < totalPageCount; i++)
        {
            PDPage page = document.getPage(i);
            int pageNumber = i + 1;
            String displayString = pageNumber + "/" + totalPageCount;

            try(PDPageContentStream cs = new PDPageContentStream(document, page, PDPageContentStream.AppendMode.APPEND, true, true))
            {
                ReportResources.TextStyle pageNumStyle = mReportResources.pageNumberStyle();
                cs.beginText();
                cs.setFont(pageNumStyle.font(), pageNumStyle.fontSize());
                cs.setNonStrokingColor(pageNumStyle.color());
                cs.newLineAtOffset(58, 20);
                cs.showText(displayString);
                cs.endText();

                if(mAddDisclaimer)
                {
                    ReportResources.TextStyle disclaimerStyle = mReportResources.disclaimerStyle();
                    String disclaimer = "This report is for 'Research Use Only (RUO)' and is not suitable "
                            + "for diagnostic or clinical applications. No rights can be derived from the contents of this report.";
                    cs.beginText();
                    cs.setFont(disclaimerStyle.font(), disclaimerStyle.fontSize());
                    cs.setNonStrokingColor(disclaimerStyle.color());
                    cs.newLineAtOffset(98, 20);
                    cs.showText(disclaimer);
                    cs.endText();
                }
            }
            catch(IOException e)
            {
                throw new RuntimeException("Failed to write footer", e);
            }
        }
    }
}
