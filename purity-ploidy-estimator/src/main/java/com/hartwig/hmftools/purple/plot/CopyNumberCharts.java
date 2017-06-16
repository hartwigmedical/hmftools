package com.hartwig.hmftools.purple.plot;

import java.awt.Color;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

class CopyNumberCharts {

    static JFreeChart cumulativePloidy(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final XYDataset dataset = CopyNumberCharts.createDataset(copyNumbers);
        JFreeChart chart = ChartFactory.createScatterPlot("Cumulative Ploidy Distribution", "BAF Weighting (CDF)",
                "Ploidy", dataset, PlotOrientation.VERTICAL, false, false, false);
        XYPlot xyPlot = (XYPlot) chart.getPlot();
        XYItemRenderer renderer = ((XYPlot) chart.getPlot()).getRenderer();
        Shape shape = new Ellipse2D.Double(0, 0, 4, 4);
        renderer.setSeriesShape(0, shape);
        renderer.setSeriesPaint(0, Color.blue);
        xyPlot.getRangeAxis().setRange(0, 10);
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

}
