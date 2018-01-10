package com.hartwig.hmftools.purple.plot;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.Comparator;
import java.util.List;
import java.util.function.ToDoubleFunction;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYBarRenderer;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.CategoryTableXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

class CopyNumberCharts {

    private static final int MAX_COPY_NUMBER_SERIES = 6;

    // JOBA: http://colorbrewer2.org/
    private static final Color COPY_NUMBER_1 = new Color(255, 26, 28);
    private static final Color COPY_NUMBER_2 = new Color(77, 175, 74);
    private static final Color COPY_NUMBER_3 = new Color(55, 126, 184);
    private static final Color COPY_NUMBER_4 = new Color(152, 78, 163);
    private static final Color COPY_NUMBER_5 = new Color(255, 127, 0);
    private static final Color COPY_NUMBER_6 = new Color(255, 255, 51);

    @NotNull
    static JFreeChart copyNumberCDF(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final XYDataset dataset = CopyNumberCharts.createDataset(copyNumbers);
        JFreeChart chart =
                ChartFactory.createScatterPlot("Copy Number CDF", "BAF Weighting (CDF)", "Ploidy", dataset, PlotOrientation.VERTICAL, false,
                        false, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        XYItemRenderer renderer = ((XYPlot) chart.getPlot()).getRenderer();
        Shape shape = new Ellipse2D.Double(0, 0, 4, 4);
        renderer.setSeriesShape(0, shape);
        renderer.setSeriesPaint(0, Color.blue);
        xyPlot.getRangeAxis().setRange(0, 10);
        return chart;
    }

    @NotNull
    static JFreeChart copyNumberPDF(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final XYDataset dataset = ploidyPDF(copyNumbers, PurpleCopyNumber::averageTumorCopyNumber, PurpleCopyNumber::bafCount);
        JFreeChart chart =
                ChartFactory.createScatterPlot("Copy Number PDF", "Ploidy", "BAF Count", dataset, PlotOrientation.VERTICAL, false, false,
                        false);

        XYBarRenderer renderer = new XYBarRenderer(0.9);
        renderer.setShadowVisible(false);
        renderer.setBarPainter(new StandardXYBarPainter());
        renderer.setSeriesPaint(0, Color.BLUE);

        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.getDomainAxis().setRange(0, 10);
        xyPlot.setRenderer(renderer);
        return chart;
    }

    @NotNull
    static JFreeChart minorAllelePDF(@NotNull final List<PurpleCopyNumber> variants) {
        final CategoryTableXYDataset dataset = minorAllele(variants);
        final JFreeChart chart =
                ChartFactory.createXYBarChart("Minor Allele Ploidy PDF", "Ploidy", false, "BAF Count", dataset, PlotOrientation.VERTICAL,
                        true, false, false);

        StackedXYBarRenderer renderer = new StackedXYBarRenderer();
        renderer.setBarPainter(new StandardXYBarPainter());
        renderer.setShadowVisible(false);

        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setRenderer(renderer);
        xyPlot.getDomainAxis().setRange(-1.1, 5.1);

        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            renderer.setSeriesPaint(i, copyNumberColor(String.valueOf(dataset.getSeriesKey(i))));
        }
        return chart;
    }

    @NotNull
    static JFreeChart somaticPloidyPDF(@NotNull final List<PurityAdjustedSomaticVariant> variants) {
        final CategoryTableXYDataset dataset = variants(variants);
        final JFreeChart chart =
                ChartFactory.createXYBarChart("Somatic Variant Ploidy PDF", "Ploidy", false, "Count", dataset, PlotOrientation.VERTICAL,
                        true, false, false);

        StackedXYBarRenderer renderer = new StackedXYBarRenderer();
        renderer.setBarPainter(new StandardXYBarPainter());
        renderer.setShadowVisible(false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        xyPlot.setRenderer(renderer);

        for (int i = 0; i < dataset.getSeriesCount(); i++) {
            renderer.setSeriesPaint(i, copyNumberColor(String.valueOf(dataset.getSeriesKey(i))));
        }
        return chart;
    }

    @NotNull
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

    @NotNull
    private static <T> XYDataset ploidyPDF(@NotNull List<T> events, @NotNull ToDoubleFunction<T> function,
            @NotNull ToDoubleFunction<T> increment) {
        final int bucketCount = 100;
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

    @NotNull
    private static CategoryTableXYDataset variants(@NotNull final List<PurityAdjustedSomaticVariant> variants) {
        int maxPloidy = 6;
        int maxBuckets = maxPloidy * 10;

        double[][] buckets = new double[MAX_COPY_NUMBER_SERIES + 1][maxBuckets];
        for (final PurityAdjustedSomaticVariant variant : variants) {
            double value = variant.adjustedVAF() * variant.adjustedCopyNumber();

            int series = (int) Math.max(0, Math.min(MAX_COPY_NUMBER_SERIES, Math.round(variant.adjustedCopyNumber())));
            int column = Math.min(maxBuckets - 1, Math.max(0, (int) Math.round(value / 0.1)));
            buckets[series][column] += 1;
        }

        CategoryTableXYDataset result = new CategoryTableXYDataset();
        for (int i = 0; i <= MAX_COPY_NUMBER_SERIES; i++) {
            String seriesName = "CN" + i + (i == MAX_COPY_NUMBER_SERIES ? "+" : "");

            for (int j = 0; j < maxBuckets; j++) {
                if (!Doubles.isZero(buckets[i][j])) {
                    result.add(j * 0.1, buckets[i][j], seriesName);
                }
            }
        }

        return result;
    }

    @VisibleForTesting
    @NotNull
    static CategoryTableXYDataset minorAllele(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        int positivePloidy = 5;
        int positiveBuckets = positivePloidy * 10;
        int negativePloidy = 1;
        int negativeBuckets = negativePloidy * 10;
        int totalBuckets = negativeBuckets + positiveBuckets;

        double[][] buckets = new double[MAX_COPY_NUMBER_SERIES + 1][totalBuckets + 1];
        for (final PurpleCopyNumber copyNumber : copyNumbers) {

            double unboundMinorAllele = (1 - copyNumber.averageActualBAF()) * copyNumber.averageTumorCopyNumber();

            int series = limit(copyNumber.averageTumorCopyNumber(), 0, MAX_COPY_NUMBER_SERIES);
            int column = limit(unboundMinorAllele / 0.1, -negativeBuckets, positiveBuckets) + negativeBuckets;

            buckets[series][column] += copyNumber.bafCount();
        }

        CategoryTableXYDataset result = new CategoryTableXYDataset();
        for (int i = 0; i <= MAX_COPY_NUMBER_SERIES; i++) {
            String seriesName = "CN" + i + (i == MAX_COPY_NUMBER_SERIES ? "+" : "");

            for (int j = 0; j <= totalBuckets; j++) {
                if (!Doubles.isZero(buckets[i][j])) {
                    result.add(-1 + (j * 0.1), buckets[i][j], seriesName);
                }
            }
        }

        return result;
    }

    private static int limit(double value, int min, int max) {
        return Math.max(min, Math.min(max, (int) Math.round(value)));
    }

    @NotNull
    private static Color copyNumberColor(final @NotNull String seriesName) {
        switch (seriesName) {
            case "CN1":
                return COPY_NUMBER_1;
            case "CN2":
                return COPY_NUMBER_2;
            case "CN3":
                return COPY_NUMBER_3;
            case "CN4":
                return COPY_NUMBER_4;
            case "CN5":
                return COPY_NUMBER_5;
            case "CN6+":
                return COPY_NUMBER_6;
        }

        return Color.BLACK;
    }
}
