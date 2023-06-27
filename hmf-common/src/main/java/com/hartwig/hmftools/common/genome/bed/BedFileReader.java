package com.hartwig.hmftools.common.genome.bed;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class BedFileReader
{
    private static final Logger LOGGER = LogManager.getLogger(BedFileReader.class);

    public static boolean loadBedFile(final String filename, final List<ChrBaseRegion> regions)
    {
        try
        {
            regions.addAll(loadBedFile(filename));
            return true;
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load BED file: {}, error: {}", filename, e.toString());
            return false;
        }
    }

    public static List<ChrBaseRegion> loadBedFile(final String filename) throws Exception
    {
        BufferedReader reader = createBufferedReader(filename);

        List<String> lines = Lists.newArrayList();
        String line = null;

        while((line = reader.readLine()) != null)
        {
            lines.add(line);
        }

        List<ChrBaseRegion> regions = loadBedFile(lines);
        LOGGER.info("loaded {} regions from BED file {}", regions.size(), filename);
        return regions;
    }

    public static List<ChrBaseRegion> loadBedFile(final List<String> lines) throws Exception
    {
        List<ChrBaseRegion> regions = Lists.newArrayList();

        for(String line : lines)
        {
            if(line.contains("Chromosome"))
                continue;

            final String[] values = line.split(TSV_DELIM, -1);

            if(values.length < 3)
            {
                throw new Exception("invalid slice BED entry: " + line);
            }

            String chromosome = values[0];
            int posStart = Integer.parseInt(values[1]) + 1;
            int posEnd = Integer.parseInt(values[2]);
            regions.add(new ChrBaseRegion(chromosome, posStart, posEnd));
        }

        return regions;
    }
}
