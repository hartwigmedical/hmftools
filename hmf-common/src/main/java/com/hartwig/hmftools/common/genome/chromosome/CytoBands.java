package com.hartwig.hmftools.common.genome.chromosome;

import static java.lang.String.format;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.FileDelimiters;

import org.apache.logging.log4j.core.util.IOUtils;

public class CytoBands
{
    private final Map<String,List<CytoBand>> mChromosomeMap;

    public CytoBands(final RefGenomeVersion refGenomeVersion)
    {
        mChromosomeMap = Maps.newHashMap();
        loadRegions(refGenomeVersion);
    }

    public @Nullable String getCytoBandName(final String chromosome, final int position)
    {
        List<CytoBand> bands = mChromosomeMap.get(chromosome);

        if(bands == null)
        {
            return null;
        }

        CytoBand band = bands.stream().filter(x -> x.containsPosition(position)).findFirst().orElse(null);
        return band != null ? band.Name : null;
    }

    private void loadRegions(final RefGenomeVersion refGenomeVersion)
    {
        String resourceFilename = resourceFilename(refGenomeVersion);

        List<String> lines = new BufferedReader(new InputStreamReader(CytoBands.class.getResourceAsStream(resourceFilename)))
                .lines().collect(Collectors.toList());

        String delim = FileDelimiters.inferFileDelimiter(resourceFilename);

        int chrIndex = 0;
        int posStartIndex = 1;
        int posEndIndex = 2;
        int nameIndex = 3;

        lines.remove(0); // remove header

        for(String line : lines)
        {
            final String[] values = line.split(delim, -1);

            String chromosome = refGenomeVersion.versionedChromosome(values[chrIndex]);
            int posStart = Integer.parseInt(values[posStartIndex]);
            int posEnd = Integer.parseInt(values[posEndIndex]);

            List<CytoBand> bands = mChromosomeMap.get(chromosome);

            if(bands == null)
            {
                bands = Lists.newArrayList();
                mChromosomeMap.put(chromosome, bands);
            }

            bands.add(new CytoBand(posStart, posEnd, values[nameIndex]));
        }
    }

    private class CytoBand extends BaseRegion
    {
        public final String Name;

        public CytoBand(final int positionStart, final int positionEnd, final String name)
        {
            super(positionStart, positionEnd);
            Name = name;
        }
    }

    // methods to get the resource file as-is
    public static String resourceAsString(final RefGenomeVersion refGenomeVersion)
    {
        try
        {
            String resourceFilename = resourceFilename(refGenomeVersion);
            return readResource(resourceFilename);
        }
        catch(IOException e)
        {
            return "";
        }
    }

    public static String resourceFilename(final RefGenomeVersion refGenomeVersion)
    {
        return format("/refgenome/cytoBands.%s.tsv", refGenomeVersion.identifier());
    }

    private static String readResource(final String resource) throws IOException
    {
        InputStream in = CytoBands.class.getResourceAsStream(resource);
        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
        return IOUtils.toString(reader);
    }
}
