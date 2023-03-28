package com.hartwig.hmftools.orange.report.components;

import java.net.URL;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

import org.jetbrains.annotations.NotNull;

public class Header {

    @NotNull
    private final PdfImageXObject orangeCircosObject;
    private final ReportResources reportResources;

    public Header(@NotNull URL orangeCircosPath, @NotNull ReportResources reportResources) {
        this.orangeCircosObject = new PdfImageXObject(ImageDataFactory.create(orangeCircosPath));
        this.reportResources = reportResources;
    }

    public void renderHeader(@NotNull PdfPage page) {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        pdfCanvas.addXObject(orangeCircosObject, 50, page.getPageSize().getHeight() - 70, 60, false);

        cv.add(new Paragraph().add(new Text("O").setFont(reportResources.fontBold())
                        .setFontSize(11)
                        .setFontColor(ReportResources.PALETTE_ORANGE_1))
                .add(new Text("R").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_2))
                .add(new Text("A").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_3))
                .add(new Text("N").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_4))
                .add(new Text("G").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_5))
                .add(new Text("E").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_6))
                .add(new Text(" Report").setFont(reportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLACK))
                .setFixedPosition(175, page.getPageSize().getHeight() - 40, 300));

        pdfCanvas.release();
    }
}
