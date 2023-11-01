package com.hartwig.hmftools.amber;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

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

        writer.write(format("chromosome\tstartPosition\tendPosition\tsnpCount\thomCount\thetCount\tfilter"));
        writer.newLine();

        for(RegionOfHomozygosity region : regions)
        {
            writer.write(format("%s\t%d\t%d\t%d\t%d\t%d\t%s",
                    region.Chromosome.toString(), region.Start, region.End, region.getSnpCount(),
                    region.NumHomozygous, region.NumHeterozygous, softFilter(region)));
            writer.newLine();
        }

        writer.close();
    }

    // we have several soft filters
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
