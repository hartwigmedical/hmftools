package com.hartwig.hmftools.common.drivercatalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class DriverCatalogFile {

    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static final String DELIMITER = "\t";
    private static final String OLD_DRIVER_CATALOG_EXTENSION = ".driver.catalog.tsv";
    private static final String DRIVER_CATALOG_EXTENSION = ".driver.catalog.somatic.tsv";
    private static final String GERMLINE_DRIVER_CATALOG_EXTENSION = ".driver.catalog.germline.tsv";

    @NotNull
    public static String generateSomaticFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = generateSomaticFilenameForWriting(basePath, sample);
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + OLD_DRIVER_CATALOG_EXTENSION;
    }

    public static String generateSomaticFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + DRIVER_CATALOG_EXTENSION;
    }

    public static String generateGermlineFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + GERMLINE_DRIVER_CATALOG_EXTENSION;
    }

    @NotNull
    public static List<DriverCatalog> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final List<DriverCatalog> catalog) throws IOException {
        Files.write(new File(filename).toPath(), toLines(catalog));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<DriverCatalog> catalog) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        catalog.stream().map(DriverCatalogFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<DriverCatalog> fromLines(@NotNull final List<String> lines) {
        return lines.stream().skip(1).map(x -> fromString(x)).collect(Collectors.toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("chromosome")
                .add("chromosomeBand")
                .add("gene")
                .add("driver")
                .add("category")
                .add("likelihoodMethod")
                .add("driverLikelihood")
                .add("NA")
                .add("missense")
                .add("nonsense")
                .add("splice")
                .add("inframe")
                .add("frameshift")
                .add("biallelic")
                .add("minCopyNumber")
                .add("maxCopyNumber")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DriverCatalog ratio) {
        return new StringJoiner(DELIMITER).add(ratio.chromosome())
                .add(ratio.chromosomeBand())
                .add(ratio.gene())
                .add(String.valueOf(ratio.driver()))
                .add(String.valueOf(ratio.category()))
                .add(String.valueOf(ratio.likelihoodMethod()))
                .add(FORMAT.format(ratio.driverLikelihood()))
                .add("0")
                .add(String.valueOf(ratio.missense()))
                .add(String.valueOf(ratio.nonsense()))
                .add(String.valueOf(ratio.splice()))
                .add(String.valueOf(ratio.inframe()))
                .add(String.valueOf(ratio.frameshift()))
                .add(String.valueOf(ratio.biallelic()))
                .add(FORMAT.format(ratio.minCopyNumber()))
                .add(FORMAT.format(ratio.maxCopyNumber()))
                .toString();
    }

    @NotNull
    private static DriverCatalog fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .chromosome(values[0])
                .chromosomeBand(values[1])
                .gene(values[2])
                .driver(DriverType.valueOf(values[3]))
                .category(DriverCategory.valueOf(values[4]))
                .likelihoodMethod(LikelihoodMethod.valueOf(values[5]))
                .driverLikelihood(Double.parseDouble(values[6]))
                .missense(Long.parseLong(values[8]))
                .nonsense(Long.parseLong(values[9]))
                .splice(Long.parseLong(values[10]))
                .inframe(Long.parseLong(values[11]))
                .frameshift(Long.parseLong(values[12]))
                .biallelic(Boolean.parseBoolean(values[13]))
                .minCopyNumber(Double.parseDouble(values[14]))
                .maxCopyNumber(Double.parseDouble(values[15]));

        return builder.build();
    }
}
