package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.renderer.DivRenderer;
import com.itextpdf.layout.renderer.DrawContext;
import com.itextpdf.layout.renderer.IRenderer;

public class InlineBarChart extends Div {

    private double value;
    private double min;
    private double max;

    public InlineBarChart(double value, double min, double max) {
        this.value = value;
        this.min = min;
        this.max = max;
    }

    public final double getValue() {
        return value;
    }

    final double getMin() {
        return min;
    }

    public final double getMax() {
        return max;
    }

    /**
     * Remap v from [inMin, inMax] to [outMin, outMax]
     */
    protected static double map(double v, double inMin, double inMax, double outMin, double outMax) {
        return (v - inMin) * (outMax - outMin) / (inMax - inMin) + outMin;
    }

    @Override
    public IRenderer getRenderer() {
        return new BarChartRenderer(this);
    }

    private class BarChartRenderer extends DivRenderer {

        BarChartRenderer(final InlineBarChart inlineBarChart) {
            super(inlineBarChart);
        }

        public void draw(final DrawContext drawContext) {
            final PdfCanvas canvas = drawContext.getCanvas();

            final Rectangle area = this.occupiedArea.getBBox();
            final float x = area.getX();
            final float y = area.getY();
            final float width = area.getWidth();
            final float height = area.getHeight();
            final float radius = height * .5f;
            final float filledWidth = (float) map(getValue(), getMin(), getMax(), height, width);

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
