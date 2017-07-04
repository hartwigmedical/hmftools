package com.hartwig.hmftools.purple.plot;

import java.awt.Color;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.SymbolAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.MatrixSeriesCollection;
import org.jfree.data.xy.NormalizedMatrixSeries;

class BubbleCharts {

    private static final double MAX_COPY_NUMBER = 5;
    private static final double MAX_VARIANT_PLOIDY = 5;
    private static final double MAX_BUBBLE_SIZE = 15;

    public static JFreeChart createCopyNumberBubblePlot(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final MatrixSeriesCollection dataset = createBubbleSet(copyNumbers);
        JFreeChart chart =
                ChartFactory.createBubbleChart("CopyNumber vs BAF", "BAF (%)", "Ploidy", dataset, PlotOrientation.VERTICAL, false, false, false);

        chart.getXYPlot().getRenderer(0).setSeriesPaint(0, Color.blue);

        XYPlot plot = (XYPlot) chart.getPlot();
        plot.getDomainAxis().setRange(46, 104);
        plot.setRangeAxis(ploidyAxis(MAX_COPY_NUMBER));

        return chart;
    }

    public static JFreeChart createVariantBubblePlot(@NotNull final List<EnrichedSomaticVariant> copyNumbers) {
        final MatrixSeriesCollection dataset = createVariantBubbleSet(copyNumbers);
        JFreeChart chart =
                ChartFactory.createBubbleChart("Variant Ploidy vs VAF ", "VAF (%)", "Variant Ploidy", dataset, PlotOrientation.VERTICAL, false, false, false);

        chart.getXYPlot().getRenderer(0).setSeriesPaint(0, new Color(140, 140, 100));

        XYPlot plot = (XYPlot) chart.getPlot();
        plot.getDomainAxis().setRange(0, 104);
        plot.setRangeAxis(ploidyAxis(MAX_VARIANT_PLOIDY));

        return chart;
    }

    private static NumberAxis ploidyAxis(double max) {
        final NumberFormat numberFormat = new DecimalFormat("##");
        List<String> labels = Lists.newArrayList();
        for (int i = 0; i <= 110; i++) {
            labels.add(String.valueOf(numberFormat.format((double) i / 10)));
        }
        final NumberAxis axis = new SymbolAxis("Ploidy", labels.toArray(new String[labels.size()]));
        axis.setTickUnit(new NumberTickUnit(10));
        axis.setRange(0, max * 10);

        return axis;
    }

    private static MatrixSeriesCollection createBubbleSet(@NotNull final List<PurpleCopyNumber> copyNumbers) {

        double maxAbsoluteValue = 0;
        int maxX = bubbleBafIndex(1);
        int maxY = bubbleCopyNumberIndex(MAX_COPY_NUMBER);

        NormalizedMatrixSeries series = new NormalizedMatrixSeries("Bubble", maxY + 1, maxX + 1);

        for (PurpleCopyNumber copyNumber : copyNumbers) {
            if (Doubles.greaterOrEqual(copyNumber.averageTumorCopyNumber(), 0) && Doubles.greaterOrEqual(copyNumber.averageActualBAF(),
                    0.4)) {
                int x = bubbleBafIndex(Math.min(1, copyNumber.averageActualBAF()));
                int y = bubbleCopyNumberIndex(Math.min(MAX_COPY_NUMBER, copyNumber.averageTumorCopyNumber()));

                double current = series.get(y, x);
                maxAbsoluteValue = Math.max(current, maxAbsoluteValue);
                series.update(y, x, current + copyNumber.bafCount());
            }
        }

        int totalWeigh = copyNumbers.stream().mapToInt(PurpleCopyNumber::bafCount).sum();
        series.setScaleFactor(MAX_BUBBLE_SIZE / maxAbsoluteValue * totalWeigh);
        return new MatrixSeriesCollection(series);
    }

    private static MatrixSeriesCollection createVariantBubbleSet(@NotNull final List<EnrichedSomaticVariant> variants) {

        double maxAbsoluteValue = 0;
        int maxX = bubbleBafIndex(1.1);
        int maxY = bubbleCopyNumberIndex(MAX_VARIANT_PLOIDY);

        NormalizedMatrixSeries series = new NormalizedMatrixSeries("Bubble", maxY + 1, maxX + 1);

        for (EnrichedSomaticVariant variant : variants) {

            int x = bubbleBafIndex(Math.max(0, Math.min(1.1, variant.adjustedVAF())));
            int y = bubbleCopyNumberIndex(Math.min(MAX_VARIANT_PLOIDY, variant.adjustedVAF() * Math.round(variant.adjustedCopyNumber())));

            double current = series.get(y, x);
            maxAbsoluteValue = Math.max(current, maxAbsoluteValue);
            series.update(y, x, current + 1);
        }

        int totalWeigh = variants.size();
        series.setScaleFactor(MAX_BUBBLE_SIZE / maxAbsoluteValue * totalWeigh);
        return new MatrixSeriesCollection(series);
    }

    static int bubbleCopyNumberIndex(double copyNumber) {
        return (int) Math.round(copyNumber * 10);
    }

    static int bubbleBafIndex(double copyNumber) {
        return (int) Math.round(copyNumber * 100);
    }

}
