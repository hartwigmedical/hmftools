package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;
import org.jetbrains.annotations.NotNull;

public final class Header {

    private static final String HMF_LOGO_PATH = "pdf/hartwig_logo.jpg";

    private PdfImageXObject hmfLogoObj;

    public Header() {

        // Attempt to load image object
        ImageData imgData = ReportResources.loadImageData(HMF_LOGO_PATH);
        if (imgData != null) {
            hmfLogoObj = new PdfImageXObject(imgData);
        }

    }

    public void renderHeader(@NotNull String chapterTitle, boolean firstPageOfChapter, @NotNull PdfPage page) {

        final PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());



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
