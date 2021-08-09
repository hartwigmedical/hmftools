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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class DriverCatalogFile
{
    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static final String DELIMITER = "\t";
    private static final String SOMATIC_DRIVER_CATALOG_EXTENSION = ".driver.catalog.somatic.tsv";
    private static final String GERMLINE_DRIVER_CATALOG_EXTENSION = ".driver.catalog.germline.tsv";

    @NotNull
    public static String generateSomaticFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + SOMATIC_DRIVER_CATALOG_EXTENSION;
    }

    @NotNull
    public static String generateGermlineFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + GERMLINE_DRIVER_CATALOG_EXTENSION;
    }

    @NotNull
    public static List<DriverCatalog> read(@NotNull final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final List<DriverCatalog> catalog) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(catalog));
    }

    @VisibleForTesting
    @NotNull
    static List<String> toLines(@NotNull final List<DriverCatalog> catalog)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        catalog.stream().map(DriverCatalogFile::toString).forEach(lines::add);
        return lines;
    }

    @VisibleForTesting
    @NotNull
    static List<DriverCatalog> fromLines(@NotNull final List<String> lines)
    {
        return lines.stream().skip(1).map(x -> fromString(x)).collect(Collectors.toList());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("chromosome")
                .add("chromosomeBand")
                .add("gene")
                .add("driver")
                .add("category")
                .add("likelihoodMethod")
                .add("driverLikelihood")
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
    private static String toString(@NotNull final DriverCatalog driverCatalog)
    {
        return new StringJoiner(DELIMITER).add(driverCatalog.chromosome())
                .add(driverCatalog.chromosomeBand())
                .add(driverCatalog.gene())
                .add(String.valueOf(driverCatalog.driver()))
                .add(String.valueOf(driverCatalog.category()))
                .add(String.valueOf(driverCatalog.likelihoodMethod()))
                .add(FORMAT.format(driverCatalog.driverLikelihood()))
                .add(String.valueOf(driverCatalog.missense()))
                .add(String.valueOf(driverCatalog.nonsense()))
                .add(String.valueOf(driverCatalog.splice()))
                .add(String.valueOf(driverCatalog.inframe()))
                .add(String.valueOf(driverCatalog.frameshift()))
                .add(String.valueOf(driverCatalog.biallelic()))
                .add(FORMAT.format(driverCatalog.minCopyNumber()))
                .add(FORMAT.format(driverCatalog.maxCopyNumber()))
                .toString();
    }

    @NotNull
    private static DriverCatalog fromString(@NotNull final String line)
    {
        // TODO: Clean up the version with 16 entries. This is the entry that contains dndsDriverLikelihood
        //      This can be cleaned up following the instructions in DEV-1924
        String[] values = line.split(DELIMITER);
        if(values.length == 16)
        {
            return ImmutableDriverCatalog.builder().chromosome(values[0])
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
                    .maxCopyNumber(Double.parseDouble(values[15])).build();
        }
        else if(values.length == 15)
        {
            return ImmutableDriverCatalog.builder().chromosome(values[0])
                    .chromosomeBand(values[1])
                    .gene(values[2])
                    .driver(DriverType.valueOf(values[3]))
                    .category(DriverCategory.valueOf(values[4]))
                    .likelihoodMethod(LikelihoodMethod.valueOf(values[5]))
                    .driverLikelihood(Double.parseDouble(values[6]))
                    .missense(Long.parseLong(values[7]))
                    .nonsense(Long.parseLong(values[8]))
                    .splice(Long.parseLong(values[9]))
                    .inframe(Long.parseLong(values[10]))
                    .frameshift(Long.parseLong(values[11]))
                    .biallelic(Boolean.parseBoolean(values[12]))
                    .minCopyNumber(Double.parseDouble(values[13]))
                    .maxCopyNumber(Double.parseDouble(values[14]))
                    .build();
        }
        else
        {
            throw new IllegalStateException("Invalid driver catalog entry found: '" + line + "'");
        }
    }
}
