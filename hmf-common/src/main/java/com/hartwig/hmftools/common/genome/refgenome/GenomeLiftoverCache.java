package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.Options;

public class GenomeLiftoverCache
{
    private final Map<String,List<CoordMapping>> mMappings;

    public static final String LIFTOVER_MAPPING_FILE = "liftover_mapping";
    public static final String LIFTOVER_MAPPING_FILE_DESC = "Liftover mapping file";

    private static final String DELIM = "\t";
    private static final String NEG_ORIENT = "-";

    public GenomeLiftoverCache()
    {
        mMappings = Maps.newHashMap();
    }

    public boolean hasMappings() { return !mMappings.isEmpty(); }

    public List<CoordMapping> getChromosomeMappings(final String chromosome) { return mMappings.get(chromosome); }

    public static void addConfig(final Options options)
    {
        options.addOption(LIFTOVER_MAPPING_FILE, true, LIFTOVER_MAPPING_FILE_DESC);
    }

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
                    LOGGER.error("invalid mapping line: {}", line);
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

            LOGGER.info("loaded {} coordinate mapping entries from file: {}", mMappings.size(), filename);
            return true;
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load coordinate mapping entries from file: {}", filename, e.toString());
            return false;
        }
    }
}
