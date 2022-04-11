package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.jetbrains.annotations.NotNull;

public final class RegionOfHomozygosityFile
{
    private static final String AMBER_EXTENSION = ".amber.homozygousregion.tsv";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull final List<RegionOfHomozygosity> regions) throws IOException
    {
        CSVFormat csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(
                        "chromosome", "startPosition", "endPosition", "snpCount", "homCount", "hetCount", "filter")
                .build();
        try (CSVPrinter csvPrinter = csvFormat.print(createBufferedWriter(filename)))
        {
            for (RegionOfHomozygosity region : regions)
            {
                csvPrinter.print(region.getChromosome().toString());
                csvPrinter.print(region.getStart());
                csvPrinter.print(region.getEnd());
                csvPrinter.print(region.getSnpCount());
                csvPrinter.print(region.getNumHomozygous());
                csvPrinter.print(region.getNumHeterozygous());
                csvPrinter.print(softFilter(region));
                csvPrinter.println();
            }
        }
    }

    // we have several soft filters
    private static String softFilter(@NotNull final RegionOfHomozygosity region)
    {
        List<String> softFilters = new ArrayList<>();

        if (region.getLength() < AmberConstants.HOMOZYGOUS_REGION_LONG_SIZE)
            softFilters.add("minLength");

        double hetRatio = ((double)region.getNumHeterozygous()) / region.getSnpCount();

        if (hetRatio > AmberConstants.HOMOZYGOUS_REGION_MAX_HET_RATIO)
            softFilters.add("maxHetProportion");

        if (!softFilters.isEmpty())
            return String.join(";", softFilters);
        return "PASS";
    }
}
