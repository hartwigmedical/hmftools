package com.hartwig.hmftools.common.drivercatalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class DriverCatalogFile {

    static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    static final String HEADER_PREFIX = "Gene";
    private static final String DELIMITER = "\t";
    private static final String DRIVER_CATALOG_EXTENSION = ".driver.catalog";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + DRIVER_CATALOG_EXTENSION;
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
        final List<DriverCatalog> result = Lists.newArrayList();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                result.add(fromString(line));
            }
        }
        return result;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("Gene")
                .add("Category")
                .add("Driver")
                .add("DriverLikelihood")
                .add("DndsLikelihood")
                .add("Missense")
                .add("Nonsense")
                .add("Splice")
                .add("Inframe")
                .add("Frameshift")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final DriverCatalog ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.gene()))
                .add(String.valueOf(ratio.category()))
                .add(String.valueOf(ratio.driver()))
                .add(FORMAT.format(ratio.driverLikelihood()))
                .add(FORMAT.format(ratio.dndsLikelihood()))
                .add(String.valueOf(ratio.missense()))
                .add(String.valueOf(ratio.nonsense()))
                .add(String.valueOf(ratio.splice()))
                .add(String.valueOf(ratio.inframe()))
                .add(String.valueOf(ratio.frameshift()))
                .toString();
    }

    @NotNull
    private static DriverCatalog fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder()
                .gene(values[0])
                .category(DriverCategory.valueOf(values[1]))
                .driver(DriverType.valueOf(values[2]))
                .driverLikelihood(Double.valueOf(values[3]))
                .dndsLikelihood(Double.valueOf(values[4]))
                .missense(Long.valueOf(values[5]))
                .nonsense(Long.valueOf(values[6]))
                .splice(Long.valueOf(values[7]))
                .inframe(Long.valueOf(values[8]))
                .frameshift(Long.valueOf(values[9]));

        return builder.build();
    }
}
