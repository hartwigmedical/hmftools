package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegion;
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

    public void write(@NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<EnrichedSomaticVariant> variants) throws IOException {

        final List<PurpleCopyNumber> filteredCopyNumber =
                copyNumbers.stream().filter(x -> !isSexChromosome(x)).filter(x -> x.bafCount() > 0).collect(Collectors.toList());

        final List<EnrichedSomaticVariant> filteredSomaticVariants =
                variants.stream().filter(x -> !isSexChromosome(x)).collect(Collectors.toList());

        copyNumberPDF(filteredCopyNumber);
        minorAllelePloidyPDF(filteredCopyNumber);
        somaticPloidyPDF(filteredSomaticVariants);
    }

    private void copyNumberCDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber_CDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.copyNumberCDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    private void copyNumberPDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.copyNumberPDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    private void somaticPloidyPDF(@NotNull final List<EnrichedSomaticVariant> variants) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".variant.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.somaticPloidyPDF(variants);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    private void minorAllelePloidyPDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".minor_allele.png";
        JFreeChart chart = CopyNumberCharts.minorAllelePDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), chart, 500, 300);
    }

    private static boolean isSexChromosome(GenomeRegion region) {
        return region.chromosome().equals("X") || region.chromosome().equals("Y");
    }

    private static boolean isSexChromosome(GenomePosition region) {
        return region.chromosome().equals("X") || region.chromosome().equals("Y");
    }

}
