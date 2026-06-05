package com.hartwig.hmftools.orange.report;

import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import java.awt.image.BufferedImage;

public final class PdfConverter
{
    public static void convertPdfToPng(final String pdfFilename, final String pngFilename) throws IOException
    {
        //        PDDocument document = PDDocument.load(new File(pdfFilename));
        //
        //        PDFRenderer pdfRenderer = new PDFRenderer(document);
        //
        //        for (int page = 0; page < document.getNumberOfPages(); ++page)
        //        {
        //            BufferedImage bim = pdfRenderer.renderImageWithDPI(page, 300); // 300 DPI = high quality
        //
        //            ImageIO.write(bim, "png", new File(pngFilename));
        //        }
        //
        //        document.close();
    }
}
