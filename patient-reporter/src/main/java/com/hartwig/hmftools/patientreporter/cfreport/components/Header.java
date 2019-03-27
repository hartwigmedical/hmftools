package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

public final class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg";
    private static PdfImageXObject hmfLogoObj = null;

    public static void addHeader(String chapterTitle, boolean firstPageOfChapter, PdfPage page) {

        final PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        // Attempt to load image object
        if (hmfLogoObj == null) {
            hmfLogoObj = new PdfImageXObject(ReportResources.loadImageData(HMF_LOGO_PATH));
        }

        // Add HMF logo image object to page
        if (hmfLogoObj != null) {
            pdfCanvas.addXObject(hmfLogoObj, 52, 772, 44, false);
        }

        // Add "Hartwig Medical OncoAct"
        cv.add(new Paragraph()
                .add(new Text("Hartwig Medical")
                        .setFont(ReportResources.getFontBold())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text(" Onco")
                        .setFont(ReportResources.getFontRegular())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_BLUE))
                .add(new Text("Act")
                        .setFont(ReportResources.getFontRegular())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_RED))
                .setFixedPosition(230, 791, 300));

        // Add chapter title
        Paragraph chapterTitleParagraph = new Paragraph(chapterTitle)
                .addStyle(ReportResources.chapterTitleStyle())
                .setFixedPosition(ReportResources.PAGE_MARGIN_LEFT, 721, 500);
        if (!firstPageOfChapter) {
            chapterTitleParagraph.add(new Text(" (Continued)").setFontSize(11));
        }

        cv.add(chapterTitleParagraph);

        pdfCanvas.release();

    }

}
