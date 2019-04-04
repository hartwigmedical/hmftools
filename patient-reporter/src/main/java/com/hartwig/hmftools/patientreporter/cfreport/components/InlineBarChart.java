package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.renderer.DivRenderer;
import com.itextpdf.layout.renderer.DrawContext;
import com.itextpdf.layout.renderer.IRenderer;

public class InlineBarChart extends Div {

    public static Scale LINEAR_SCALE = (v) -> { return v; };
    public static Scale LOG10_SCALE = (v) -> { return Math.log10(v); };

    private double value;
    private double min;
    private double max;

    private Scale scale = LINEAR_SCALE; // Default to linear scale

    public InlineBarChart(double value, double min, double max) {
        this.value = value;
        this.min = min;
        this.max = max;
    }

    public void setScale(Scale scale) {
        this.scale = scale;
    }

    public final double getValue() {
        return value;
    }

    public final double getScaledValue() {
        return getScaledValue(getValue());
    }

    public final double getMin() {
        return min;
    }

    public final double getScaledMin() {
        return getScaledValue(getMin());
    }

    public final double getMax() {
        return max;
    }

    public final double getScaledMax() {
        return getScaledValue(getMax());
    }

    protected final double getScaledValue(double value) {
        return scale.scale(value);
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
            final float filledWidth = (float) MathUtil.map(getScaledValue(), getScaledMin(), getScaledMax(), height, width);

            // Background
            canvas.setFillColor(ReportResources.PALETTE_LIGHT_GREY);
            canvas.roundRectangle(x, y, width, height, radius);
            canvas.fill();

            // Fill
            if (getValue() > getMin()) {
                canvas.setFillColor(ReportResources.PALETTE_BLUE);
                canvas.roundRectangle(x, y, filledWidth, height, radius);
                canvas.fill();
            }

            super.draw(drawContext);

        }

    }

    public interface Scale {
        double scale(final double value);
    }

}
