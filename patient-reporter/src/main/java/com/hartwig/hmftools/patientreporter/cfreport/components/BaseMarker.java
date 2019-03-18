package com.hartwig.hmftools.patientreporter.cfreport.components;


import com.hartwig.hmftools.patientreporter.cfreport.ReportConfiguration;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfWriter;
import org.jetbrains.annotations.NotNull;

import java.awt.*;
import java.util.Random;

/**
 * Class to draw one or more rounded rectangles (DNA Base reference)
 */
public class BaseMarker {

    private static final float WIDTH = 12.3f;
    private static final float HEIGHT = 1.9f;
    private static final float RADIUS = HEIGHT * .5f;

    private static final float LINE_WIDTH = .25f;

    /**
     * Draw line of base markers at the bottom of the page
     * @param fullWidth
     * @param pdfWriter
     */
    public static void drawFooterLine(boolean fullWidth, @NotNull PdfWriter pdfWriter) {

        // Positioning in pt
        final float startX = 156f;
        final float y = 22;
        final float deltaX = 87;

        final float markerCount = fullWidth ? 5 : 3;
        final float redProbability = .2f;

        PdfContentByte pb = pdfWriter.getDirectContent();
        pb.saveState();

        Random r = new Random();
        for (int i = 0; i < markerCount; i++) {
            float x = startX + i * deltaX;
            Color c = (r.nextFloat() < redProbability) ? ReportConfiguration.PALETTE_RED : ReportConfiguration.PALETTE_CYAN;
            BaseMarker.drawMarker(x, y, c, false, pb);
        }

        pb.restoreState();


    }

    public static void drawSidepanelMarkers(boolean fullHeight, @NotNull PdfWriter pdfWriter) {

        // Positioning in pt
        final float startX = 455;
        final float startY = fullHeight  ? 22 : 775;
        final float deltaX = 35;
        final float deltaY = 42;

        final float redProbability = .05f;
        final float filledProbability = .15f;

        final int rowCount = fullHeight ? 20 : 2;

        PdfContentByte pb = pdfWriter.getDirectContent();
        pb.saveState();

        Random r = new Random();
        for (int row = 0; row < rowCount; row++) {
            for (int col = 0; col < 4; col++) {

                float x = startX + col * deltaX;
                float y = startY + row * deltaY;
                Color c = (r.nextFloat() < redProbability) ? ReportConfiguration.PALETTE_RED : ReportConfiguration.PALETTE_CYAN;
                boolean filled = (r.nextFloat() < filledProbability);

                BaseMarker.drawMarker(x, y, c, filled, pb);

            }

        }

        pb.restoreState();

    }

    private static void drawMarker(float x, float y, @NotNull Color color, boolean filled, @NotNull PdfContentByte contentByte) {

        contentByte.roundRectangle(x, y, WIDTH, HEIGHT, RADIUS);
        contentByte.setLineWidth(LINE_WIDTH);

        if (filled) {

            // Apply both fill and stroke
            contentByte.setColorStroke(color);
            contentByte.setColorFill(color);
            contentByte.fillStroke();

        } else {

            // Apply only stroke
            contentByte.setColorStroke(color);
            contentByte.stroke();

        }

    }

}
