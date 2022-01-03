package com.hartwig.hmftools.common.drivercatalog;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.HOM_DEL_DISRUPTION;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.checkConvertType;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

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

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("chromosome")
                .add("chromosomeBand")
                .add("gene")
                .add("transcript")
                .add("isCanonical")
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
                .add(driverCatalog.transcript())
                .add(String.valueOf(driverCatalog.isCanonical()))
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
    public static List<DriverCatalog> read(@NotNull final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @VisibleForTesting
    @NotNull
    static List<DriverCatalog> fromLines(@NotNull final List<String> lines)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            drivers.add(ImmutableDriverCatalog.builder()
                    .chromosome(values[fieldsIndexMap.get("chromosome")])
                    .chromosomeBand(values[fieldsIndexMap.get("chromosomeBand")])
                    .gene(values[fieldsIndexMap.get("gene")])
                    .transcript(fieldsIndexMap.containsKey("transcript") ?
                            values[fieldsIndexMap.get("transcript")] : "")
                    .isCanonical(fieldsIndexMap.containsKey("isCanonical") ?
                            Boolean.parseBoolean(values[fieldsIndexMap.get("isCanonical")]) : true)
                    .driver(checkConvertType(values[fieldsIndexMap.get("driver")]))
                    .category(DriverCategory.valueOf(values[fieldsIndexMap.get("category")]))
                    .likelihoodMethod(LikelihoodMethod.valueOf(values[fieldsIndexMap.get("likelihoodMethod")]))
                    .driverLikelihood(Double.parseDouble(values[fieldsIndexMap.get("driverLikelihood")]))
                    .missense(Integer.parseInt(values[fieldsIndexMap.get("missense")]))
                    .nonsense(Integer.parseInt(values[fieldsIndexMap.get("nonsense")]))
                    .splice(Integer.parseInt(values[fieldsIndexMap.get("splice")]))
                    .inframe(Integer.parseInt(values[fieldsIndexMap.get("inframe")]))
                    .frameshift(Integer.parseInt(values[fieldsIndexMap.get("frameshift")]))
                    .biallelic(Boolean.parseBoolean(values[fieldsIndexMap.get("biallelic")]))
                    .minCopyNumber(Double.parseDouble(values[fieldsIndexMap.get("minCopyNumber")]))
                    .maxCopyNumber(Double.parseDouble(values[fieldsIndexMap.get("maxCopyNumber")]))
                    .build());
        }

        return drivers;
    }
}
