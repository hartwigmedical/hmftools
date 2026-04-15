package com.hartwig.hmftools.orange.report.components;

import java.io.IOException;

import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.font.PDFont;

public class SidePanel
{
    private static final float ROW_SPACING = 35;
    private static final float VALUE_TEXT_Y_OFFSET = 18;
    private static final float MAX_WIDTH = 120;

    private static final float RECTANGLE_WIDTH = 170;
    private static final float RECTANGLE_HEIGHT = 60;

    private final String mSampleId;
    private final ReportResources mReportResources;

    public SidePanel(final String sampleId, final ReportResources reportResources)
    {
        mSampleId = sampleId;
        mReportResources = reportResources;
    }

    public void renderSidePanel(final PDPage page, final PDDocument document)
    {
        try(PDPageContentStream cs = new PDPageContentStream(document, page, PDPageContentStream.AppendMode.APPEND, true, true))
        {
            PDRectangle pageSize = page.getMediaBox();

            // Draw orange rectangle in top-right corner
            float rectX = pageSize.getWidth() - RECTANGLE_WIDTH;
            float rectY = pageSize.getHeight() - RECTANGLE_HEIGHT;
            cs.setNonStrokingColor(ReportResources.PALETTE_ORANGE);
            cs.addRect(rectX, rectY, RECTANGLE_WIDTH, RECTANGLE_HEIGHT);
            cs.fill();

            // Draw "SAMPLE" label
            float xPos = pageSize.getWidth() - RECTANGLE_WIDTH + 15;
            float yPos = (pageSize.getHeight() + 15) - ROW_SPACING;

            ReportResources.TextStyle labelStyle = mReportResources.sidePanelLabelStyle();
            cs.beginText();
            cs.setFont(labelStyle.font(), labelStyle.fontSize());
            cs.setNonStrokingColor(labelStyle.color());
            cs.newLineAtOffset(xPos, yPos);
            cs.showText("SAMPLE");
            cs.endText();

            // Draw sample ID value
            float valueFontSize = maxPointSizeForWidth(mReportResources.fontBold(), 11, 6, mSampleId, MAX_WIDTH);
            yPos -= VALUE_TEXT_Y_OFFSET;

            cs.beginText();
            cs.setFont(mReportResources.fontBold(), valueFontSize);
            cs.setNonStrokingColor(ReportResources.PALETTE_WHITE);
            cs.newLineAtOffset(xPos, yPos);
            cs.showText(mSampleId);
            cs.endText();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to render side panel", e);
        }
    }

    private static float maxPointSizeForWidth(
            final PDFont font, float initialFontSize, float minFontSize, final String text, float maxWidth)
    {
        float fontIncrement = 0.1F;
        float fontSize = initialFontSize;

        try
        {
            float width = font.getStringWidth(text) / 1000 * fontSize;
            while(width > maxWidth && fontSize > minFontSize)
            {
                fontSize -= fontIncrement;
                width = font.getStringWidth(text) / 1000 * fontSize;
            }
        }
        catch(IOException e)
        {
            // Fall back to initial font size
        }

        return fontSize;
    }
}
