package com.hartwig.hmftools.orange.e2e;

import java.awt.Rectangle;
import java.io.IOException;

import com.google.common.base.Preconditions;

import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.text.PDFTextStripperByArea;

public class PageTextRipper
{
    private final PDPage page;
    private final int pageHeight;
    private final int pageWidth;

    public PageTextRipper(final PDPage page)
    {
        this.page = page;
        pageHeight = (int) page.getMediaBox().getHeight();
        pageWidth = (int) page.getMediaBox().getWidth();
    }

    public String[] getLinesInRectangle(PositionInPage topLeft, PositionInPage bottomRight) throws IOException
    {
        Preconditions.checkArgument(topLeft.fractionAcrossPage < bottomRight.fractionAcrossPage);
        Preconditions.checkArgument(topLeft.fractionDownPage < bottomRight.fractionDownPage);
        int topLeftX = (int) (topLeft.fractionAcrossPage * pageWidth);
        int topLeftY = (int) (topLeft.fractionDownPage * pageHeight);
        int width = (int) ((bottomRight.fractionAcrossPage - topLeft.fractionAcrossPage) * pageWidth);
        int height = (int) ((bottomRight.fractionDownPage - topLeft.fractionDownPage) * pageHeight);
        Rectangle region = new Rectangle(topLeftX, topLeftY, width, height);

        PDFTextStripperByArea byArea = new PDFTextStripperByArea();
        byArea.setSortByPosition(true);
        byArea.addRegion("region", region);
        byArea.extractRegions(page);
        return byArea.getTextForRegion("region").split("\n");
    }
}
