package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.kernel.pdf.canvas.PdfCanvas;
import com.itextpdf.layout.Canvas;
import com.itextpdf.layout.element.Div;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.TextAlignment;
import com.itextpdf.layout.renderer.DivRenderer;
import com.itextpdf.layout.renderer.DrawContext;
import com.itextpdf.layout.renderer.IRenderer;
import org.jetbrains.annotations.NotNull;

import java.text.NumberFormat;

public class BarChart extends InlineBarChart {

    private final static float HEIGHT = 45;

    private String lowLabel;
    private String highLabel;
    private String[] tickLabels = {};

    private boolean overshootEnabled = false;
    private String overshootLabel = "";

    private Indicator indicator = null;


    public BarChart(float value, float min, float max, @NotNull String lowLabel, @NotNull String highLabel, @NotNull String[] tickLabels) {
        super(value, min, max);
        super.setHeight(HEIGHT);
        this.lowLabel = lowLabel;
        this.highLabel = highLabel;
        this.tickLabels = tickLabels;
    }

    public void setIndicator(float value, @NotNull String name) {
        if (value >= getMin() && value <= getMax()) {
            indicator = new Indicator(name, value);
        } else {
            System.err.println("Indicator value outside bounds");
        }
    }

    public void enableDefaultRangeOvershoot(@NotNull String overshootLabel) {
        overshootEnabled = true;
        this.overshootLabel = overshootLabel;
    }

    @Override
    public Div setHeight(float height) {
        System.err.println("Cannot set height of BarChart, has fixed height");
        return super.setHeight(HEIGHT);
    }

    @Override
    public IRenderer getRenderer() {
        return new BarChartRenderer(this);
    }

    public final static String[] createTickMarkLabels(float min, float max, float increment, @NotNull NumberFormat format) {

        int valueCount = 1 + (int) ((max - min) / increment);

        String[] tickMarkLabels = new String[valueCount];
        for (int i = 0; i < valueCount; i++) {
            tickMarkLabels[i] = format.format(min + i * increment);
        }

        return tickMarkLabels;

    }

    private class BarChartRenderer extends DivRenderer {

        private final BarChart barChart;

        private final float BAR_OUTLINE_HEIGHT = 9f;
        private final float BAR_INSET = 2f;

        public BarChartRenderer(final BarChart barChart) {
            super(barChart);
            this.barChart = barChart;
        }

        public void draw(final DrawContext drawContext) {
            super.draw(drawContext);

            final PdfCanvas canvas = drawContext.getCanvas();
            final Rectangle boundingBox = this.occupiedArea.getBBox();
            final Canvas cv = new Canvas(canvas, canvas.getDocument(), boundingBox);

            final float valueBarWidth = overshootEnabled ? boundingBox.getWidth() - 40 : boundingBox.getWidth();

            final float clampedValue = Math.max(getMin(), Math.min(getMax(), getValue()));

            final Rectangle barOutlineRect = new Rectangle(boundingBox.getX(), boundingBox.getTop() - 39, valueBarWidth, BAR_OUTLINE_HEIGHT);

            final float outerBarRadius = barOutlineRect.getHeight() * .5f;

            // Overshoot indicator
            if (overshootEnabled) {

                final float barOverlap = 2;

                // Overshoot outline
                final Rectangle overshootOutlineRect = new Rectangle(barOutlineRect.getX(), barOutlineRect.getY(), boundingBox.getWidth(), barOutlineRect.getHeight());
                canvas.setStrokeColor(ReportResources.PALETTE_MID_BLUE);
                canvas.setLineWidth(.25f);
                canvas.setLineDash(3f, 2f);
                canvas.roundRectangle(overshootOutlineRect.getX(), overshootOutlineRect.getY(), overshootOutlineRect.getWidth(), overshootOutlineRect.getHeight(), outerBarRadius);
                canvas.stroke();
                canvas.setLineDash(1f);


                // Filled bar
                canvas.setFillColor(ReportResources.PALETTE_BLUE);
                if (getValue() > getMax()) {
                    final Rectangle barRect = new Rectangle(overshootOutlineRect.getX() + BAR_INSET, overshootOutlineRect.getY() +  BAR_INSET, overshootOutlineRect.getWidth() - 2 * BAR_INSET, overshootOutlineRect.getHeight() - 2 * BAR_INSET);
                    canvas.roundRectangle(barRect.getX(), barRect.getY(), barRect.getWidth(), barRect.getHeight(), barRect.getHeight() * .5);
                    canvas.fill();
                }

                // White mask
                canvas.setFillColor(ReportResources.PALETTE_WHITE);
                canvas.circle(barOutlineRect.getRight() - outerBarRadius, barOutlineRect.getY() + outerBarRadius, outerBarRadius + barOverlap);
                canvas.fill();

                // Overshoot label
                cv.showTextAligned(new Paragraph(overshootLabel)
                        .addStyle(ReportResources.subTextStyle().setFontSize(6)), overshootOutlineRect.getRight() - 4, barOutlineRect.getBottom() - 21f, TextAlignment.RIGHT);

            }

            // Outer border
            canvas.setStrokeColor(ReportResources.PALETTE_BLUE);
            canvas.setFillColor(ReportResources.PALETTE_WHITE);
            canvas.setLineWidth(.25f);
            canvas.roundRectangle(barOutlineRect.getX(), barOutlineRect.getY(), barOutlineRect.getWidth(), barOutlineRect.getHeight(), outerBarRadius);
            canvas.fillStroke();

            // Filled bar
            final Rectangle barRect = new Rectangle(barOutlineRect.getX() + BAR_INSET, barOutlineRect.getY() +  BAR_INSET, barOutlineRect.getWidth() - 2 * BAR_INSET, barOutlineRect.getHeight() - 2 * BAR_INSET);
            canvas.setFillColor(ReportResources.PALETTE_BLUE);
            canvas.roundRectangle(barRect.getX(), barRect.getY(), map(clampedValue, getMin(), getMax(), barRect.getHeight(), barRect.getWidth()), barRect.getHeight(), barRect.getHeight() * .5);
            canvas.fill();

            // Add top labels
            cv.showTextAligned(new Paragraph(lowLabel)
                    .addStyle(ReportResources.smallBodyHeadingStyle()), boundingBox.getLeft(), boundingBox.getTop() - 25, TextAlignment.LEFT);
            cv.showTextAligned(new Paragraph(highLabel)
                    .addStyle(ReportResources.smallBodyHeadingStyle()), boundingBox.getRight(), boundingBox.getTop() - 25, TextAlignment.RIGHT);

            // Add tick marks
            final float tickMarkDelta = (barRect.getWidth() - barRect.getHeight()) / ((float) tickLabels.length - 1);
            for (int i = 0; i < tickLabels.length ; i++) {

                float x = barRect.getLeft() + barRect.getHeight() * 0.5f + i * tickMarkDelta;

                canvas.moveTo(x, barOutlineRect.getBottom() - 4.1f);
                canvas.lineTo(x, barOutlineRect.getBottom() - 9.4f);
                canvas.setLineWidth(.25f);
                canvas.setStrokeColor(ReportResources.PALETTE_BLACK);
                canvas.stroke();

                cv.showTextAligned(new Paragraph(tickLabels[i])
                        .addStyle(ReportResources.subTextStyle().setFontSize(6)), x, barOutlineRect.getBottom() - 21f, TextAlignment.CENTER);

            }

            // Add boundary indicator
            if (indicator != null) {

                float x = map(indicator.value, getMin(), getMax(), barRect.getLeft(), barRect.getRight());

                canvas.moveTo(x, barOutlineRect.getTop() + 9.5f);
                canvas.lineTo(x, barOutlineRect.getBottom() - 10.5f);
                canvas.setLineWidth(1f);
                canvas.setStrokeColor(ReportResources.PALETTE_PINK);
                canvas.stroke();

                cv.showTextAligned(new Paragraph("\u2192 " + indicator.name.toUpperCase())
                        .addStyle(ReportResources.subTextBoldStyle().setFontSize(6)).setFontColor(ReportResources.PALETTE_PINK), x + 4.5f, barOutlineRect.getTop() + 3.5f, TextAlignment.LEFT);

            }

        }

    }

    private class Indicator {

        public String name;
        public float value;

        public Indicator(@NotNull String name, float value) {
            this.name = name;
            this.value = value;
        }

    }

}
