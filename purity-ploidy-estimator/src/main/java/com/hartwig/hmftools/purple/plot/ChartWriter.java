package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

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

    public void copyNumberCDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber_CDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.copyNumberCDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    public void copyNumberPDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber_PDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.copyNumberPDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    public void somaticPloidy(@NotNull final List<EnrichedSomaticVariant> variants)  throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".variant_PDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.somaticPloidyPDF(variants);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

}
