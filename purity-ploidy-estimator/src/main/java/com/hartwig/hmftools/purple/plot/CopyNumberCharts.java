package com.hartwig.hmftools.purple.plot;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.Comparator;
import java.util.List;
import java.util.function.ToDoubleFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

class CopyNumberCharts {

    static JFreeChart copyNumberCDF(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final XYDataset dataset = CopyNumberCharts.createDataset(copyNumbers);
        JFreeChart chart = ChartFactory.createScatterPlot("Copy Number CDF", "BAF Weighting (CDF)", "Ploidy", dataset, PlotOrientation.VERTICAL, false, false,
                false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        XYItemRenderer renderer = ((XYPlot) chart.getPlot()).getRenderer();
        Shape shape = new Ellipse2D.Double(0, 0, 4, 4);
        renderer.setSeriesShape(0, shape);
        renderer.setSeriesPaint(0, Color.blue);
        xyPlot.getRangeAxis().setRange(0, 10);
        return chart;
    }

    static JFreeChart copyNumberPDF(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final XYDataset dataset = ploidyPDF(100, copyNumbers, PurpleCopyNumber::averageTumorCopyNumber, PurpleCopyNumber::bafCount);
        JFreeChart chart = ChartFactory.createScatterPlot("Copy Number PDF", "Ploidy", "BAF Count", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        XYItemRenderer renderer = new XYLineAndShapeRenderer();
        Shape shape = new Ellipse2D.Double(-0.5, -0.5, 0, 0);
        renderer.setSeriesShape(0, shape);
        xyPlot.setRenderer(renderer);
        renderer.setSeriesPaint(0, Color.blue);
        return chart;
    }

    static JFreeChart somaticPloidyPDF(@NotNull final List<EnrichedSomaticVariant> variants) {
        final XYDataset dataset = createVariantsPDF(variants);
        JFreeChart chart = ChartFactory.createScatterPlot("Somatic Variant Ploidy PDF", "Ploidy", "Count", dataset, PlotOrientation.VERTICAL, false, false,
                false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        XYItemRenderer renderer = new XYLineAndShapeRenderer();
        Shape shape = new Ellipse2D.Double(-0.5, -0.5, 0, 0);
        renderer.setSeriesShape(0, shape);
        xyPlot.setRenderer(renderer);
        renderer.setSeriesPaint(0, new Color(140, 140, 100));
        return chart;
    }

    private static XYDataset createDataset(@NotNull final List<PurpleCopyNumber> copyNumbers) {

        final List<PurpleCopyNumber> sortedCopyNumbers = Lists.newArrayList(copyNumbers);
        sortedCopyNumbers.sort(Comparator.comparingDouble(PurpleCopyNumber::averageTumorCopyNumber));

        int totalCount = copyNumbers.stream().mapToInt(PurpleCopyNumber::bafCount).sum();

        XYSeries series = new XYSeries("CDF");
        int cumulativeBAFCount = 0;
        for (PurpleCopyNumber sortedCopyNumber : sortedCopyNumbers) {
            cumulativeBAFCount += sortedCopyNumber.bafCount();
            series.add((double) cumulativeBAFCount / totalCount * 100, sortedCopyNumber.averageTumorCopyNumber());
        }

        return new XYSeriesCollection(series);
    }

    private static XYDataset createVariantsPDF(@NotNull final List<EnrichedSomaticVariant> variants) {
        return ploidyPDF(55, variants, (EnrichedSomaticVariant x) -> x.adjustedVAF() * x.adjustedCopyNumber(), x -> 1);
    }

    private static <T> XYDataset ploidyPDF(int bucketCount, List<T> events, ToDoubleFunction<T> function, ToDoubleFunction<T> increment) {

        double[] buckets = new double[bucketCount];
        for (T event : events) {
            double value = function.applyAsDouble(event);
            int index = Math.min(bucketCount - 1, Math.max(0, (int) Math.round(value / 0.1)));
            buckets[index] = buckets[index] + increment.applyAsDouble(event);
        }

        XYSeries series = new XYSeries("PDF");
        for (int i = 0; i < bucketCount; i++) {
            series.add(i * 0.1, buckets[i]);
        }

        return new XYSeriesCollection(series);
    }

}
