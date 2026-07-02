package com.hartwig.hmftools.orange.report.components;

import static com.hartwig.hmftools.orange.report.ReportResources.HEADER_ORANGE_HEIGHT;

import java.awt.Color;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;

public class Header
{
    private final PDImageXObject mCircosImage;
    private final ReportResources mReportResources;
    private final boolean mAddDisclaimer;

    public Header(final URL orangeCircosPath, final ReportResources reportResources, boolean addDisclaimer,
            final PDDocument document)
    {
        mReportResources = reportResources;
        mAddDisclaimer = addDisclaimer;
        try(InputStream imageStream = orangeCircosPath.openStream())
        {
            mCircosImage = PDImageXObject.createFromByteArray(document, imageStream.readAllBytes(), "orange_circos.png");
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to load header image", e);
        }
    }

    public void renderHeader(final PDPage page, final PDDocument document)
    {
        try(PDPageContentStream cs = new PDPageContentStream(document, page, PDPageContentStream.AppendMode.APPEND, true, true))
        {
            PDRectangle pageSize = page.getMediaBox();
            float pageHeight = pageSize.getHeight();

            // Draw the circos image (reuses the single cached PDImageXObject across all pages)
            float imgScale = (float) HEADER_ORANGE_HEIGHT / mCircosImage.getHeight();
            float imgWidth = mCircosImage.getWidth() * imgScale;
            float imgHeight = HEADER_ORANGE_HEIGHT;
            float orangeImageVerticalPosition = pageHeight - HEADER_ORANGE_HEIGHT - 10;
            cs.drawImage(mCircosImage, 50, orangeImageVerticalPosition, imgWidth, imgHeight);

            // Draw the "ORANGE Report" title text, each letter in a different orange shade
            float fontSize = 11;
            float left = mAddDisclaimer ? 150 : 180;
            float textY = pageHeight - 40;

            Color[] letterColors = {
                    ReportResources.PALETTE_ORANGE_1, ReportResources.PALETTE_ORANGE_2,
                    ReportResources.PALETTE_ORANGE_3, ReportResources.PALETTE_ORANGE_4,
                    ReportResources.PALETTE_ORANGE_5, ReportResources.PALETTE_ORANGE_6
            };
            String[] letters = { "O", "R", "A", "N", "G", "E" };

            cs.beginText();
            cs.setFont(mReportResources.fontBold(), fontSize);
            cs.newLineAtOffset(left, textY);

            for(int i = 0; i < letters.length; i++)
            {
                cs.setNonStrokingColor(letterColors[i]);
                cs.showText(letters[i]);
            }

            cs.setNonStrokingColor(ReportResources.PALETTE_BLACK);
            cs.showText(" Report");

            if(mAddDisclaimer)
            {
                cs.setFont(mReportResources.fontBold(), 9);
                cs.showText(" (Research Use Only)");
            }

            cs.endText();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to render header", e);
        }
    }
}
