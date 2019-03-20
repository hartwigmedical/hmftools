package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.renderer.DivRenderer;
import com.itextpdf.layout.renderer.DrawContext;
import com.itextpdf.layout.renderer.IRenderer;

public class InlineBarChart extends Div {

    private float fillPercentage;

    public InlineBarChart(float value, float min, float max) {
        fillPercentage = normalize(value, min, max);
    }

    private static final float normalize(float v, float inMin, float inMax) {
        return (v - inMin) * (1f - 0f) / (inMax - inMin);
    }

    @Override
    public IRenderer getRenderer() {
        return new BarChartRenderer(this);
    }

    private class BarChartRenderer extends DivRenderer {

        private final InlineBarChart barChart;

        public BarChartRenderer(final InlineBarChart inlineBarChart) {
            super(inlineBarChart);
            this.barChart = inlineBarChart;
        }

        public void draw(final DrawContext drawContext) {
            final PdfCanvas canvas = drawContext.getCanvas();

            final Rectangle area = this.occupiedArea.getBBox();
            final float x = area.getX();
            final float y = area.getY();
            final float width = area.getWidth();
            final float height = area.getHeight();
            final float radius = height * .5f;
            final float filledWidth = barChart.fillPercentage * width;

            // Background
            canvas.setFillColor(ReportResources.PALETTE_LIGHT_GREY);
            canvas.roundRectangle(x, y, width, height, radius);
            canvas.fill();

            // Fill
            canvas.setFillColor(ReportResources.PALETTE_BLUE);
            canvas.roundRectangle(x, y, filledWidth, height, radius);
            canvas.fill();

            super.draw(drawContext);

        }

    }

}
