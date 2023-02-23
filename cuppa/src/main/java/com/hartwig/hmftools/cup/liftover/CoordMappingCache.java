package com.hartwig.hmftools.cup.liftover;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class CoordMappingCache
{
    private final Map<String,List<CoordMapping>> mMappings;

    private static final String DELIM = "\t";
    private static final String NEG_ORIENT = "-";

    public CoordMappingCache()
    {
        mMappings = Maps.newHashMap();
    }

    public boolean hasMappings() { return !mMappings.isEmpty(); }

    public List<CoordMapping> getChromosomeMappings(final String chromosome) { return mMappings.get(chromosome); }

    /*
    public void addMapping(final String chromosome, int sourceStart, int sourceEnd, int destStart, int destEnd, boolean reverse)
    {
        mMappings.add(new CoordMapping(chromosome, sourceStart, sourceEnd, destStart, destEnd, reverse));
    }
    */

    public boolean loadFile(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            lines.remove(0);

            List<CoordMapping> chrMappings = null;
            String currentChromosome = "";

            for(String line : lines)
            {
                String[] values = line.split(DELIM, -1);

                if(values.length < 8)
                {
                    CUP_LOGGER.error("invalid mapping line: {}", line);
                    return false;
                }

                String chromosome = values[0];

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    chrMappings = Lists.newArrayList();
                    mMappings.put(chromosome, chrMappings);
                }

                // note +1 for start positions since source file is in BED style
                chrMappings.add(new CoordMapping(
                        values[0], Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2]),
                        Integer.parseInt(values[5]) + 1, Integer.parseInt(values[6]), values[7].equals(NEG_ORIENT)));
            }

            CUP_LOGGER.info("loaded {} coordinate mapping entries from file: {}", mMappings.size(), filename);
            return true;
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("failed to load coordinate mapping entries from file: {}", filename, e.toString());
            return false;
        }
    }

}
