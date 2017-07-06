package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.title.TextTitle;

public class ChartWriter {

    @NotNull
    private final String sample;
    @NotNull
    private final String outputDirectory;

    public ChartWriter(@NotNull final String sample, @NotNull final String outputDirectory) {
        this.sample = sample;
        this.outputDirectory = outputDirectory;
    }

    public void write(@NotNull final FittedPurity purity, @NotNull FittedPurityScore score,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<EnrichedSomaticVariant> variants) throws IOException {

        final List<PurpleCopyNumber> filteredCopyNumber =
                copyNumbers.stream().filter(x -> !isSexChromosome(x)).filter(x -> x.bafCount() > 0).collect(Collectors.toList());

        final List<EnrichedSomaticVariant> filteredSomaticVariants =
                variants.stream().filter(x -> !isSexChromosome(x)).collect(Collectors.toList());

        final String subtitle = subtitle(sample, purity, score);

        copyNumberPDF(subtitle, filteredCopyNumber);
        minorAllelePloidyPDF(subtitle, filteredCopyNumber);
        somaticPloidyPDF(subtitle, filteredSomaticVariants);
    }

    static String subtitle(@NotNull final String sample, @NotNull final FittedPurity purity, @NotNull final FittedPurityScore score) {
        return String.format("%s PUR:%.0f%% (%.0f%%-%.0f%%) PLE:%.2f (%.2f-%.2f)",
                sample,
                purity.purity() * 100,
                score.minPurity() * 100,
                score.maxPurity() * 100,
                purity.ploidy(),
                score.minPloidy(),
                score.maxPloidy());
    }

    static String subtitle(@NotNull final String sample, final double purity, final double ploidy) {
        return String.format("%s PUR:%.0f%% PLE:%.2f", sample, purity * 100, ploidy);
    }

    private void copyNumberCDF(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber_CDF.png";
        JFreeChart cumulativePloidyChart = CopyNumberCharts.copyNumberCDF(copyNumbers);
        ChartUtilities.saveChartAsPNG(new File(fileName), cumulativePloidyChart, 500, 300);
    }

    private void copyNumberPDF(@NotNull final String subtitle, @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".copyNumber.png";
        JFreeChart chart = CopyNumberCharts.copyNumberPDF(copyNumbers);
        chart.addSubtitle(new TextTitle(subtitle));
        ChartUtilities.saveChartAsPNG(new File(fileName), chart, 500, 300);
    }

    private void somaticPloidyPDF(@NotNull final String subtitle, @NotNull final List<EnrichedSomaticVariant> variants) throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".variant.png";
        JFreeChart chart = CopyNumberCharts.somaticPloidyPDF(variants);
        chart.addSubtitle(new TextTitle(subtitle));
        ChartUtilities.saveChartAsPNG(new File(fileName), chart, 500, 300);
    }

    private void minorAllelePloidyPDF(@NotNull final String subtitle, @NotNull final List<PurpleCopyNumber> copyNumbers)
            throws IOException {
        String fileName = outputDirectory + File.separator + sample + ".minor_allele.png";
        JFreeChart chart = CopyNumberCharts.minorAllelePDF(copyNumbers);
        chart.addSubtitle(new TextTitle(subtitle));
        ChartUtilities.saveChartAsPNG(new File(fileName), chart, 500, 300);
    }

    private static boolean isSexChromosome(GenomeRegion region) {
        return region.chromosome().equals("X") || region.chromosome().equals("Y");
    }

    private static boolean isSexChromosome(GenomePosition region) {
        return region.chromosome().equals("X") || region.chromosome().equals("Y");
    }

}
