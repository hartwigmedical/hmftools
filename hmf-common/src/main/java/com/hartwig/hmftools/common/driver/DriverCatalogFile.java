package com.hartwig.hmftools.common.driver;

import static com.hartwig.hmftools.common.driver.DriverType.checkConvertType;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.ReportedStatus;

import org.jetbrains.annotations.NotNull;

public final class DriverCatalogFile
{
    static final DecimalFormat FORMAT = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.ENGLISH));

    private static final String SOMATIC_EXTENSION = ".purple.driver.catalog.somatic.tsv";
    private static final String GERMLINE_EXTENSION = ".purple.driver.catalog.germline.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isSomatic)
    {
        // backwards compatible for reading
        return basePath + File.separator + sample + (isSomatic ? SOMATIC_EXTENSION : GERMLINE_EXTENSION);
    }

    public static String generateFilenameForWriting(final String basePath, final String sample, boolean isSomatic)
    {
        return basePath + File.separator + sample + (isSomatic ? SOMATIC_EXTENSION : GERMLINE_EXTENSION);
    }

    public static String generateSomaticFilename(final String basePath, final String sample)
    {
        return generateFilename(basePath, sample, true);
    }

    public static String generateGermlineFilename(final String basePath, final String sample)
    {
        return generateFilename(basePath, sample, false);
    }

    public static void write(final String filename, final List<DriverCatalog> catalog) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(catalog));
    }

    @VisibleForTesting
    static List<String> toLines(final List<DriverCatalog> catalog)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        catalog.stream().map(DriverCatalogFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("chromosome")
                .add("chromosomeBand")
                .add("gene")
                .add("transcript")
                .add("isCanonical")
                .add("driver")
                .add("category")
                .add("likelihoodMethod")
                .add("reportedStatus")
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

    private static String toString(final DriverCatalog driverCatalog)
    {
        return new StringJoiner(TSV_DELIM)
                .add(driverCatalog.chromosome())
                .add(driverCatalog.chromosomeBand())
                .add(driverCatalog.gene())
                .add(driverCatalog.transcript())
                .add(String.valueOf(driverCatalog.isCanonical()))
                .add(String.valueOf(driverCatalog.driver()))
                .add(String.valueOf(driverCatalog.category()))
                .add(String.valueOf(driverCatalog.likelihoodMethod()))
                .add(String.valueOf(driverCatalog.reportedStatus()))
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

    public static List<DriverCatalog> read(final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @VisibleForTesting
    static List<DriverCatalog> fromLines(final List<String> lines)
    {
        List<DriverCatalog> drivers = Lists.newArrayList();

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        Integer reportedIndex = fieldsIndexMap.get("reportedStatus");

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            ReportedStatus reportedStatus = reportedIndex != null ? ReportedStatus.valueOf(values[reportedIndex]) : ReportedStatus.REPORTED;

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
                    .reportedStatus(reportedStatus)
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
