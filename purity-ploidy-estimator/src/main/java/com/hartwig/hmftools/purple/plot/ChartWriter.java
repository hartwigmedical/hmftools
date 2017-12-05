package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.title.TextTitle;

public class ChartWriter {

    @NotNull
    private final String sample;
    @NotNull
    private final String outputDirectory;

    public ChartWriter(@NotNull final String sample, @NotNull final String outputDirectory) throws IOException {
        this.sample = sample;
        this.outputDirectory = outputDirectory;

        final File output = new File(outputDirectory);
        if (!output.exists() && !output.mkdirs()) {
            throw new IOException("Unable to create plot directory " + outputDirectory);
        }
    }

    public void write(@NotNull final FittedPurity purity, @NotNull FittedPurityScore score,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<PurityAdjustedSomaticVariant> variants)
            throws IOException {
        final List<PurpleCopyNumber> filteredCopyNumber =
                copyNumbers.stream().filter(ChartWriter::isAutosome).filter(x -> x.bafCount() > 0).collect(Collectors.toList());

        final List<PurityAdjustedSomaticVariant> filteredSomaticVariants =
                variants.stream().filter(ChartWriter::isAutosome).collect(Collectors.toList());

        final String subtitle = subtitle(sample, purity, score);

        copyNumberPDF(subtitle, filteredCopyNumber);
        minorAllelePloidyPDF(subtitle, filteredCopyNumber);
        somaticPloidyPDF(subtitle, filteredSomaticVariants);
    }

    @NotNull
    private static String subtitle(@NotNull final String sample, @NotNull final FittedPurity purity,
            @NotNull final FittedPurityScore score) {
        return String.format("%s PUR:%.0f%% (%.0f%%-%.0f%%) PLE:%.2f (%.2f-%.2f)", sample, purity.purity() * 100, score.minPurity() * 100,
                score.maxPurity() * 100, purity.ploidy(), score.minPloidy(), score.maxPloidy());
    }

    @NotNull
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

    private void somaticPloidyPDF(@NotNull final String subtitle, @NotNull final List<PurityAdjustedSomaticVariant> variants)
            throws IOException {
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

    private static boolean isAutosome(@NotNull GenomeRegion region) {
        return HumanChromosome.contains(region.chromosome()) && HumanChromosome.valueOf(region).isAutosome();
    }

    private static boolean isAutosome(@NotNull GenomePosition position) {
        return HumanChromosome.contains(position.chromosome()) && HumanChromosome.valueOf(position).isAutosome();
    }
}
