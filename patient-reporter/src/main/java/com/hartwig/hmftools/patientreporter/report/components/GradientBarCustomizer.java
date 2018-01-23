package com.hartwig.hmftools.patientreporter.report.components;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.GradientPaint;
import java.util.Optional;

import org.apache.commons.lang3.tuple.Pair;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.CategoryPointerAnnotation;
import org.jfree.chart.annotations.CategoryTextAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.ui.GradientPaintTransformType;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.StandardGradientPaintTransformer;
import org.jfree.ui.TextAnchor;

import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.definition.chart.DRIChartCustomizer;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
abstract class GradientBarCustomizer implements DRIChartCustomizer {
    @NotNull
    abstract GradientPaint gradientPaint();

    @NotNull
    abstract String startText();

    @NotNull
    abstract String endText();

    abstract int value();

    @NotNull
    @Value.Default
    Optional<Pair<Integer, String>> marker() {
        return Optional.empty();
    }

    @Override
    public void customize(JFreeChart chart, ReportParameters reportParameters) {
        customizeRenderer((BarRenderer) chart.getCategoryPlot().getRenderer());
        chart.getLegend().setVisible(false);
        customizePlotOrientationAndAxis(chart.getCategoryPlot());
        marker().ifPresent(pair -> addValueMarker(chart.getCategoryPlot(), pair.getKey(), pair.getValue()));
        addValuePointer(chart.getCategoryPlot(), value());
        addStartEndLabels(chart.getCategoryPlot(), startText(), endText());
    }

    private static void addValuePointer(@NotNull final CategoryPlot categoryPlot, final int position) {
        final Object category = categoryPlot.getCategories().get(0);
        if (category instanceof Comparable) {
            final Comparable categoryKey = (Comparable) category;
            final CategoryPointerAnnotation upArrow = new CategoryPointerAnnotation("", categoryKey, position, Math.PI / 2);
            final CategoryPointerAnnotation downArrow = new CategoryPointerAnnotation("", categoryKey, position, 3 * Math.PI / 2);
            upArrow.setArrowWidth(10);
            downArrow.setArrowWidth(10);
            categoryPlot.addAnnotation(upArrow);
            categoryPlot.addAnnotation(downArrow);
        }
    }

    private static void addValueMarker(@NotNull final CategoryPlot categoryPlot, final int position, @NotNull final String text) {
        final Marker marker = new ValueMarker(position);
        final float[] dash = { 3 };
        marker.setStroke(new BasicStroke(1.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10, dash, 0));
        marker.setPaint(new Color(49, 49, 49));
        marker.setLabel(text);
        marker.setLabelAnchor(RectangleAnchor.TOP_LEFT);
        marker.setLabelTextAnchor(TextAnchor.TOP_RIGHT);
        categoryPlot.addRangeMarker(marker);
    }

    private void customizeRenderer(@NotNull final BarRenderer renderer) {
        renderer.setShadowPaint(Color.LIGHT_GRAY);
        renderer.setShadowVisible(true);
        renderer.setSeriesPaint(0, gradientPaint());
        renderer.setGradientPaintTransformer(new StandardGradientPaintTransformer(GradientPaintTransformType.HORIZONTAL));
    }

    private static void customizePlotOrientationAndAxis(@NotNull final CategoryPlot categoryPlot) {
        categoryPlot.setOrientation(PlotOrientation.HORIZONTAL);
        categoryPlot.getDomainAxis().setVisible(false);
        categoryPlot.getRangeAxis().setVisible(false);
        final ValueAxis rangeAxis = categoryPlot.getRangeAxis();
        rangeAxis.setLowerBound(0);
        rangeAxis.setUpperBound(100);
    }

    private void addStartEndLabels(@NotNull final CategoryPlot categoryPlot, @NotNull final String startText,
            @NotNull final String endText) {
        final int startPosition = (int) Math.ceil(startText.length() / 2.0) + 2;
        final int endPosition = 100 - (int) Math.ceil(endText.length() / 2.0) - 2;
        addLabel(categoryPlot, startText, startPosition, determineContrastingColor(gradientPaint().getColor1()));
        addLabel(categoryPlot, endText, endPosition, determineContrastingColor(gradientPaint().getColor2()));
    }

    private static void addLabel(@NotNull final CategoryPlot categoryPlot, @NotNull final String labelText, final int position,
            @NotNull final Color color) {
        final Object category = categoryPlot.getCategories().get(0);
        if (category instanceof Comparable) {
            final Comparable categoryKey = (Comparable) category;
            final CategoryTextAnnotation label = new CategoryTextAnnotation(labelText, categoryKey, position);
            label.setPaint(color);
            categoryPlot.addAnnotation(label);
        }
    }

    @NotNull
    private Color determineContrastingColor(@NotNull final Color color) {
        final double colorLuminance =
                0.2126 * relativeLuminance(color.getRed()) + 0.7152 * relativeLuminance(color.getGreen()) + 0.0722 * relativeLuminance(
                        color.getBlue());
        //MIVO: w3c recommended threshold for relative luminance = 0.179
        if (colorLuminance > 0.220) {
            return Color.BLACK;
        } else {
            return Color.WHITE;
        }
    }

    private double relativeLuminance(final int value) {
        final double percent = (double) value / 255.0;
        if (percent <= 0.03928) {
            return percent / 12.92;
        } else {
            return Math.pow((percent + 0.055) / 1.055, 2.4);
        }
    }
}
