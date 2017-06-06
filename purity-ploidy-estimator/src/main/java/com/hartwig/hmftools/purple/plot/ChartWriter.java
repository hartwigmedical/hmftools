package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;

public class ChartWriter {

    @NotNull
    private final String sample;
    @NotNull
    private final String outputDirectory;

    public ChartWriter(@NotNull final String sample, @NotNull final String outputDirectory) {
        this.sample = sample;
        this.outputDirectory = outputDirectory;
    }

    public void cumulativePloidy(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumberCDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.cumulativePloidy(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

}
