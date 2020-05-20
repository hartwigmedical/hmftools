package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.renderer.DivRenderer;
import com.itextpdf.layout.renderer.DrawContext;
import com.itextpdf.layout.renderer.IRenderer;

import org.jetbrains.annotations.NotNull;

public class InlineBarChart extends Div {

    public static final Scale LOG10_SCALE = Math::log10;
    private static final Scale LINEAR_SCALE = (v) -> v;

    private final double value;
    private final double min;
    private final double max;

    private boolean isEnabled = true;

    @NotNull
    private Scale scale = LINEAR_SCALE;

    public InlineBarChart(double value, double min, double max) {
        this.value = value;
        this.min = min;
        this.max = max;
    }

    public void scale(@NotNull Scale scale) {
        this.scale = scale;
    }

    final double value() {
        return value;
    }

    final double scaledValue() {
        return scaledValue(value());
    }

    public final double min() {
        return min;
    }

    final double scaledMin() {
        return scaledValue(min());
    }

    public final double max() {
        return max;
    }

    final double scaledMax() {
        return scaledValue(max());
    }

    final double scaledValue(double value) {
        return scale.scale(value);
    }

    public void enabled(boolean enabled) {
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

        BarChartRenderer(@NotNull InlineBarChart inlineBarChart) {
            super(inlineBarChart);
        }

        public void draw(@NotNull DrawContext drawContext) {
            final PdfCanvas canvas = drawContext.getCanvas();

            final Rectangle area = this.occupiedArea.getBBox();
            final float x = area.getX();
            final float y = area.getY();
            final float width = area.getWidth();
            final float height = area.getHeight();
            final float radius = height * .5f;
            final float filledWidth = (float) MathUtil.map(scaledValue(), scaledMin(), scaledMax(), height, width);

            canvas.setFillColor(ReportResources.PALETTE_LIGHT_GREY);
            canvas.roundRectangle(x, y, width, height, radius);
            canvas.fill();

            if (isEnabled() && value() > min()) {
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
