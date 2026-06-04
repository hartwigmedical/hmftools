package com.hartwig.hmftools.amber;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public final class RegionOfHomozygosityFile
{
    private static final String AMBER_EXTENSION = ".amber.homozygousregion.tsv";

    public static String generateFilename(@NotNull final String basePath, final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    public static void write(final String filename, final List<RegionOfHomozygosity> regions) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename);

        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("chromosome").add("startPosition").add("endPosition").add("snpCount").add("homCount").add("hetCount").add("filter");
        writer.write(sj.toString());
        writer.newLine();

        for(RegionOfHomozygosity region : regions)
        {
            sj = new StringJoiner(TSV_DELIM);
            sj.add(region.Chromosome.toString());
            sj.add(String.valueOf(region.Start));
            sj.add(String.valueOf(region.End));
            sj.add(String.valueOf(region.getSnpCount()));
            sj.add(String.valueOf(region.NumHomozygous));
            sj.add(String.valueOf(region.NumHeterozygous));
            sj.add(softFilter(region));

            writer.write(sj.toString());
            writer.newLine();
        }

        writer.close();
    }

    private static String softFilter(final RegionOfHomozygosity region)
    {
        List<String> softFilters = new ArrayList<>();

        if(region.getLength() < AmberConstants.HOMOZYGOUS_REGION_LONG_SIZE)
            softFilters.add("minLength");

        double hetRatio = ((double)region.NumHeterozygous) / region.getSnpCount();

        if(hetRatio > AmberConstants.HOMOZYGOUS_REGION_MAX_HET_RATIO)
            softFilters.add("maxHetProportion");

        if(!softFilters.isEmpty())
            return String.join(";", softFilters);

        return "PASS";
    }
}
