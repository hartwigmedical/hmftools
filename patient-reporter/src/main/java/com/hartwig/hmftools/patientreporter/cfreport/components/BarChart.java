package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.kernel.colors.Color;
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

    private boolean underShootEnabled = false;
    private String undershootLabel = "";

    private boolean overshootEnabled = false;
    private String overshootLabel = "";

    private Indicator[] tickMarks = {};
    private Indicator threshold = null;

    public BarChart(double value, double min, double max, @NotNull String lowLabel, @NotNull String highLabel) {
        super(value, min, max);
        super.setHeight(HEIGHT);
        this.lowLabel = lowLabel;
        this.highLabel = highLabel;
    }

    public void setTickMarks(@NotNull Indicator... tickMarks) {
        this.tickMarks = tickMarks;
    }

    public void setTickMarks(double[] values, @NotNull NumberFormat format) {
        Indicator[] tickMarks = new Indicator[values.length];
        for (int i = 0; i < values.length; i++) {
            double value = values[i];
            tickMarks[i] = new Indicator(format.format(value), value);
        }
        setTickMarks(tickMarks);
    }

    public void setTickMarks(double min, double max, double increment, @NotNull NumberFormat format) {

        int valueCount = 1 + (int) ((max - min) / increment);
        double[] values = new double[valueCount];
        for (int i = 0; i < valueCount; i++) {
            values[i] = min + i * increment;
        }

        setTickMarks(values, format);

    }

    public void setIndicator(double value, @NotNull String name) {
        if (value >= getMin() && value <= getMax()) {
            threshold = new Indicator(name, value);
        } else {
            System.err.println("Indicator value outside bounds");
        }
    }

    public void enableOvershoot(@NotNull String overshootLabel) {
        overshootEnabled = true;
        this.overshootLabel = overshootLabel;
    }

    public void enableUndershoot(@NotNull String undershootLabel) {
        underShootEnabled = true;
        this.undershootLabel = undershootLabel;
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

    private class BarChartRenderer extends DivRenderer {

        private final float BAR_OUTLINE_HEIGHT = 9f;
        private final float BAR_INSET = 2f;

        private final float OVER_UNDER_SHOOT_WIDTH = 25;
        private final float OVER_UNDERSHOOT_OVERLAP = 2;
        private final float OVER_UNDER_SHOOT_LABEL_OFFSET = 4;

        BarChartRenderer(final BarChart barChart) {
            super(barChart);
        }

        public void draw(final DrawContext drawContext) {
            super.draw(drawContext);

            // Get canvas properties
            final PdfCanvas canvas = drawContext.getCanvas();
            final Rectangle boundingBox = this.occupiedArea.getBBox();
            final Canvas cv = new Canvas(canvas, canvas.getDocument(), boundingBox);

            final float barY = boundingBox.getTop() - 39;
            final float tickY = barY - 21;

            final boolean hasUnderShoot = isEnabled() && underShootEnabled && getValue() < getMin();
            final boolean hasOverShoot = isEnabled() && overshootEnabled && getValue() > getMax();

            final Color outlineColor = isEnabled()
                    ? ReportResources.PALETTE_MID_BLUE
                    : ReportResources.PALETTE_LIGHT_GREY;
            final Color labelColor = isEnabled()
                    ? ReportResources.PALETTE_BLACK
                    : ReportResources.PALETTE_LIGHT_GREY;


            // Draw top labels
            cv.showTextAligned(new Paragraph(lowLabel)
                    .addStyle(ReportResources.smallBodyHeadingStyle().setFontColor(labelColor)), boundingBox.getLeft(), boundingBox.getTop() - 25, TextAlignment.LEFT);
            cv.showTextAligned(new Paragraph(highLabel)
                    .addStyle(ReportResources.smallBodyHeadingStyle().setFontColor(labelColor)), boundingBox.getRight(), boundingBox.getTop() - 25, TextAlignment.RIGHT);

            // Draw undershoot
            if (hasUnderShoot) {

                final float fillValue = isEnabled()
                        ? (getValue() > getMin() ? 1 : 0.1f)
                        : 0f;

                // Draw bar
                final Rectangle outerBB = new Rectangle(boundingBox.getLeft(), barY, OVER_UNDER_SHOOT_WIDTH + BAR_OUTLINE_HEIGHT, BAR_OUTLINE_HEIGHT);
                final Rectangle innerBB = getInnerRectangle(outerBB);
                drawRoundedRect(outerBB, innerBB, fillValue, outlineColor, ReportResources.PALETTE_BLUE, true, canvas);

                // Draw mask
                final float outerRadius = getHeightRadius(outerBB);
                canvas.circle(outerBB.getRight() - outerRadius, outerBB.getY() + outerRadius, outerRadius + OVER_UNDERSHOOT_OVERLAP);
                canvas.setFillColor(ReportResources.PALETTE_WHITE);
                canvas.fill();

                // Label
                cv.showTextAligned(new Paragraph(undershootLabel)
                        .addStyle(ReportResources.subTextStyle()
                                .setFontSize(6)
                                .setFontColor(labelColor)), outerBB.getLeft() + OVER_UNDER_SHOOT_LABEL_OFFSET, tickY, TextAlignment.LEFT);

            }

            // Draw overshoot
            if (hasOverShoot) {

                final float fillValue = (isEnabled() && getValue() > getMax()) ? 1 : 0;

                // Draw bar
                final Rectangle outerBB = new Rectangle(boundingBox.getRight() - OVER_UNDER_SHOOT_WIDTH - BAR_OUTLINE_HEIGHT, barY, OVER_UNDER_SHOOT_WIDTH + BAR_OUTLINE_HEIGHT, BAR_OUTLINE_HEIGHT);
                final Rectangle innerBB = getInnerRectangle(outerBB);
                drawRoundedRect(outerBB, innerBB, fillValue, outlineColor, ReportResources.PALETTE_BLUE, true, canvas);

                // Draw mask
                final float outerRadius = getHeightRadius(outerBB);
                canvas.setFillColor(ReportResources.PALETTE_WHITE);
                canvas.circle(outerBB.getLeft() + outerRadius, outerBB.getY() + outerRadius, outerRadius + OVER_UNDERSHOOT_OVERLAP);
                canvas.fill();

                // Label
                cv.showTextAligned(new Paragraph(overshootLabel)
                        .addStyle(ReportResources.subTextStyle()
                                .setFontSize(6)
                                .setFontColor(labelColor)), outerBB.getRight() - OVER_UNDER_SHOOT_LABEL_OFFSET, tickY, TextAlignment.RIGHT);

            }

            // Draw main bar
            final float fillValue = isEnabled()
                    ? (float) MathUtil.mapClamped(getScaledValue(), getScaledMin(), getScaledMax(), 0, 1)
                    : 0;

            final Rectangle mainOuterBB = new Rectangle(
                    boundingBox.getLeft() + (hasUnderShoot ? OVER_UNDER_SHOOT_WIDTH : 0),
                    barY,
                    boundingBox.getWidth()
                            - (hasUnderShoot ? OVER_UNDER_SHOOT_WIDTH : 0)
                            - (hasOverShoot ? OVER_UNDER_SHOOT_WIDTH : 0),
                    BAR_OUTLINE_HEIGHT);
            final Rectangle mainInnerBB = getInnerRectangle(mainOuterBB);
            drawRoundedRect(mainOuterBB, mainInnerBB, fillValue, outlineColor, ReportResources.PALETTE_BLUE, false, canvas);

            // Draw tick marks
            for (Indicator tickMark : tickMarks) {

                float x = (float) MathUtil.map(getScaledValue(tickMark.value), getScaledMin(), getScaledMax(), mainInnerBB.getLeft(), mainInnerBB.getRight());

                canvas.moveTo(x, tickY + 17);
                canvas.lineTo(x, tickY + 10);
                canvas.setLineWidth(.25f);
                canvas.setStrokeColor(labelColor);
                canvas.stroke();

                TextAlignment alignment = !tickMark.equals(tickMarks[tickMarks.length - 1]) ? TextAlignment.CENTER : TextAlignment.RIGHT;

                cv.showTextAligned(new Paragraph(tickMark.name)
                        .addStyle(ReportResources.subTextStyle()
                                .setFontSize(6)
                                .setFontColor(labelColor)), x, tickY, alignment);

            }

            // Draw threshold indicator
            if (isEnabled() && threshold != null) {

                float x = (float) MathUtil.map(getScaledValue(threshold.value), getScaledMin(), getScaledMax(), mainInnerBB.getLeft(), mainInnerBB.getRight());

                canvas.moveTo(x, mainOuterBB.getTop() + 9.5f);
                canvas.lineTo(x, mainOuterBB.getBottom() - 10.5f);
                canvas.setLineWidth(1f);
                canvas.setStrokeColor(ReportResources.PALETTE_PINK);
                canvas.stroke();

                cv.showTextAligned(new Paragraph("\u2192 " + threshold.name.toUpperCase())
                        .addStyle(ReportResources.subTextBoldStyle().setFontSize(6)).setFontColor(ReportResources.PALETTE_PINK), x + 4.5f, mainOuterBB.getTop() + 3.5f, TextAlignment.LEFT);

            }

        }

        private void drawRoundedRect(@NotNull Rectangle outerBoundingBox, @NotNull Rectangle innerBoundingBox, float filledPercentage, @NotNull Color outlineColor, @NotNull Color fillColor, boolean dashedOutline, @NotNull PdfCanvas canvas) {

            filledPercentage = (float) MathUtil.clamp(filledPercentage, 0, 1);
            final Rectangle outerBarBoundingBox = outerBoundingBox;

            // Outline
            canvas.setStrokeColor(outlineColor);
            canvas.setFillColor(ReportResources.PALETTE_WHITE);
            canvas.setLineWidth(.25f);
            if  (dashedOutline) {
                canvas.setLineDash(3f, 2f);
            }
            canvas.roundRectangle(outerBarBoundingBox.getX(),
                    outerBarBoundingBox.getY(),
                    outerBarBoundingBox.getWidth(),
                    outerBarBoundingBox.getHeight(),
                    getHeightRadius(outerBarBoundingBox));
            canvas.fillStroke();
            canvas.setLineDash(1f); // Reset dash

            // Filled bar
            if (filledPercentage > 0) {
                final float innerBarRadius = getHeightRadius(innerBoundingBox);
                canvas.setFillColor(fillColor);
                canvas.roundRectangle(
                        innerBoundingBox.getX(),
                        innerBoundingBox.getY(),
                        MathUtil.map(filledPercentage, 0, 1, innerBarRadius, innerBoundingBox.getWidth()),
                        innerBoundingBox.getHeight(),
                        innerBarRadius);
                canvas.fill();
            }

        }

        @NotNull
        private Rectangle getInnerRectangle(@NotNull Rectangle outerBoundingBox) {
            return new Rectangle(
                    outerBoundingBox.getX() + BAR_INSET,
                    outerBoundingBox.getY() + BAR_INSET,
                    outerBoundingBox.getWidth() - BAR_INSET * 2,
                    outerBoundingBox.getHeight() - BAR_INSET * 2);
        }

        private float getHeightRadius(@NotNull Rectangle rect) {
            return rect.getHeight() * 0.5f;
        }

    }

    public class Indicator {

        public String name;
        public double value;

        public Indicator(@NotNull String name, double value) {
            this.name = name;
            this.value = value;
        }

    }

}
