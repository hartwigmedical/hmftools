package com.hartwig.hmftools.orange.report.component;

import java.net.MalformedURLException;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.kernel.pdf.xobject.PdfFormXObject;
import com.itextpdf.kernel.pdf.xobject.PdfImageXObject;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class Header {

    private static final Logger LOGGER = LogManager.getLogger(Header.class);

    @Nullable
    private final PdfImageXObject orangeCircosObj;

    public Header(@NotNull String orangeCircosPath) {
        ImageData orangeImage = null;
        try {
            orangeImage = ImageDataFactory.create(orangeCircosPath);
        } catch (MalformedURLException e) {
            LOGGER.warn("Could not orange image from {}", orangeCircosPath);
        }
        orangeCircosObj = orangeImage != null ? new PdfImageXObject(orangeImage) : null;
    }

    public void renderHeader(@NotNull PdfPage page) {
        PdfCanvas pdfCanvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        Canvas cv = new Canvas(pdfCanvas, page.getDocument(), page.getPageSize());

        if (orangeCircosObj != null) {
            pdfCanvas.addXObject(orangeCircosObj, 52, 772, 60, false);
        }

        cv.add(new Paragraph().add(new Text("O").setFont(ReportResources.fontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_ORANGE_1))
                .add(new Text("R").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_2))
                .add(new Text("A").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_3))
                .add(new Text("N").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_4))
                .add(new Text("G").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_5))
                .add(new Text("E").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_ORANGE_6))
                .add(new Text(" Report").setFont(ReportResources.fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLACK))
                .setFixedPosition(200, 800, 300));

        pdfCanvas.addXObject(new PdfFormXObject(new Rectangle(0, 0, 200, 0)), ReportResources.PAGE_MARGIN_LEFT, 200);

        pdfCanvas.release();
    }
}
