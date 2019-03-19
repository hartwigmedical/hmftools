package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;

public final class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg"; // "src/main/resources/pdf/hartwig_logo.jpg";
    private static PdfImageXObject hmfLogoObj = null;

    public static void addHeader(String chapterTitle, PdfPage page) {

        final PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());

        // Attempt to load image object
        if (hmfLogoObj == null) {
            hmfLogoObj = new PdfImageXObject(ReportResources.loadImageData(HMF_LOGO_PATH));
        }

        // Add image object to page
        if (hmfLogoObj != null) {
            pdfCanvas.addXObject(hmfLogoObj, 52, 772, 44, false);
        }

        // @TODO: add "Hartwig Medical OncoAct"

        // Add chapter title
        Paragraph chapterTitleParagraph = new Paragraph(chapterTitle)
                .addStyle(ReportResources.chapterTitleStyle())
                .setFixedPosition(ReportResources.PAGE_MARGIN_LEFT, 721, 500);

        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());
        cv.add(chapterTitleParagraph);

        pdfCanvas.release();

    }

}
