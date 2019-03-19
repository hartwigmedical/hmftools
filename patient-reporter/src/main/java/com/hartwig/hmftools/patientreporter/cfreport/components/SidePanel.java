package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.PageEventHandler;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.PdfPage;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import org.jetbrains.annotations.NotNull;

public final class SidePanel {

    public static void addSidePanel(PdfPage page, boolean fullHeight, boolean fullContent) {

        // Draw background and markers
        final PdfCanvas canvas = new PdfCanvas(page.getLastContentStream(), page.getResources(), page.getDocument());
        final Rectangle pageSize = page.getPageSize();

        drawBackgroundRect(fullHeight, canvas, pageSize);
        BaseMarker.drawMarkerGrid(4, (fullHeight ? 20 : 2),455, 35, fullHeight ? 22 : 775, 42, .05f, .15f, canvas);
        canvas.release();

        // @TODO Add sidepanel content from SampleReport
        if (fullHeight && fullContent) {

        } else {

        }

    }

    /**
     * Draw background rectangle, either full height or only top
     *
     * @param fullHeight
     * @param canvas
     * @param pageSize
     */
    private static void drawBackgroundRect(boolean fullHeight, @NotNull PdfCanvas canvas, @NotNull Rectangle pageSize) {

        final float RECTANGLE_WIDTH = 170;           // Width of the blue rectangle in pt
        final float RECTANGLE_HEIGHT_SHORT = 110;    // Height of the blue rectangle in pt when not full page height

        canvas.rectangle(pageSize.getWidth(), pageSize.getHeight(), -RECTANGLE_WIDTH, fullHeight ? -pageSize.getHeight() : -RECTANGLE_HEIGHT_SHORT);
        canvas.setFillColor(ReportResources.PALETTE_BLUE);
        canvas.fill();

    }

 }
