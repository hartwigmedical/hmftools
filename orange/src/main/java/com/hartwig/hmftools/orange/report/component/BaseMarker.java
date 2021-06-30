package com.hartwig.hmftools.orange.report.component;

import java.util.Random;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;

import org.jetbrains.annotations.NotNull;

final class BaseMarker {

    private BaseMarker() {
    }

    static void renderMarkerGrid(float xCount, float yCount, float xStart, float xSpacing, float yStart, float ySpacing,
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

    private static void renderMarker(float x, float y, @NotNull DeviceRgb color, boolean filled, @NotNull PdfCanvas canvas) {
        float height = 1.9f;
        canvas.roundRectangle(x, y, 12.3f, height, height * .5f);
        canvas.setLineWidth(.25f);

        if (filled) {
            canvas.setFillColor(color);
            canvas.setStrokeColor(color);
            canvas.fillStroke();
        } else {
            canvas.setStrokeColor(color);
            canvas.stroke();
        }
    }
}
