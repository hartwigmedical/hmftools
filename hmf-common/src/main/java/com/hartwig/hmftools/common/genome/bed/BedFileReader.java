package com.hartwig.hmftools.common.genome.bed;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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

    public static Map<Chromosome,List<BaseRegion>> loadBedFileChrMap(final String filename)
    {
        final Map<Chromosome,List<BaseRegion>> chrRegionMap = Maps.newHashMap();

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String line = "";

            List<BaseRegion> chrRegions = null;
            Chromosome currentChr = null;

            while((line = fileReader.readLine()) != null)
            {
                if(line.contains("Chromosome"))
                    continue;

                final String[] values = line.split(TSV_DELIM, -1);

                Chromosome chromosome = HumanChromosome.fromString(values[0]);
                int posStart = Integer.parseInt(values[1]) + 1; // as per convention
                int posEnd = Integer.parseInt(values[2]);

                if(currentChr != chromosome)
                {
                    currentChr = chromosome;
                    chrRegions = Lists.newArrayList();
                    chrRegionMap.put(chromosome, chrRegions);
                }

                chrRegions.add(new BaseRegion(posStart, posEnd));
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load panel BED file({}): {}", filename, e.toString());
            return null;
        }

        return chrRegionMap;
    }

}
