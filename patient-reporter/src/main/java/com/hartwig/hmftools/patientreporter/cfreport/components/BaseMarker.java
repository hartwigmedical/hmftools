package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;

import org.jetbrains.annotations.NotNull;

import java.util.Random;

/**
 * Class to draw a rounded rectangle (as DNA Base reference) for both the footer and the
 * side panel
 */
class BaseMarker {

    /**
     * Draw single marker at x,y with given color. PdfCanvas is *not* released after drawing
     */
    private static void renderMarker(float x, float y, @NotNull DeviceRgb color, boolean filled, @NotNull PdfCanvas canvas) {

        final float height = 1.9f;
        canvas.roundRectangle(x, y, 12.3f, height, height * .5f);
        canvas.setLineWidth(.25f);

        if (filled) {

            // Apply both fill and stroke
            canvas.setFillColor(color);
            canvas.setStrokeColor(color);
            canvas.fillStroke();

        } else {

            // Apply only stroke
            canvas.setStrokeColor(color);
            canvas.stroke();

        }

    }

    /**
     * Draw pattern grid of markers. PdfCanvas is *not* released after drawing
     *
     * @param xCount            number of columns
     * @param yCount            number of rows
     * @param xStart            horizontal start of markers in pt
     * @param xSpacing          horizontal offset between markers in pt
     * @param yStart            vertical start of markers in pt
     * @param ySpacing          vertical offset between markers in pt
     * @param redProbability    probability of a red marker [0.0, 1.0]; higher value is more probable
     * @param filledProbability probability of a filled marker [0.0, 1.0]; higher value is more probable
     * @param canvas            canvas to draw on
     */
    public static void renderMarkerGrid(float xCount, float yCount, float xStart, float xSpacing, float yStart, float ySpacing,
            float redProbability, float filledProbability, @NotNull PdfCanvas canvas) {

        Random r = new Random();
        for (int row = 0; row < yCount; row++) {

            for (int col = 0; col < xCount; col++) {

                float x = xStart + col * xSpacing;
                float y = yStart + row * ySpacing;
                DeviceRgb color = (r.nextFloat() < redProbability) ? ReportResources.PALETTE_RED : ReportResources.PALETTE_CYAN;
                boolean filled = (r.nextFloat() < filledProbability);

                renderMarker(x, y, color, filled, canvas);

            }

        }

    }

}


