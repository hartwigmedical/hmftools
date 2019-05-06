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

    private static Scale LINEAR_SCALE = (v) -> v;
    public static Scale LOG10_SCALE = Math::log10;

    private double value;
    private double min;
    private double max;

    private boolean isEnabled = true;

    private Scale scale = LINEAR_SCALE; // Default to linear scale

    public InlineBarChart(double value, double min, double max) {
        this.value = value;
        this.min = min;
        this.max = max;
    }

    public void setScale(Scale scale) {
        this.scale = scale;
    }

    final double getValue() {
        return value;
    }

    final double getScaledValue() {
        return getScaledValue(getValue());
    }

    public final double getMin() {
        return min;
    }

    final double getScaledMin() {
        return getScaledValue(getMin());
    }

    public final double getMax() {
        return max;
    }

    final double getScaledMax() {
        return getScaledValue(getMax());
    }

    final double getScaledValue(double value) {
        return scale.scale(value);
    }

    public void setEnabled(boolean enabled) {
        isEnabled = enabled;
    }

    boolean isEnabled() {
        return isEnabled;
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
            if (isEnabled() && getValue() > getMin()) {
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
