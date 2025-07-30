package com.hartwig.hmftools.common.region;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.function.Function;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class BedLine
{
    public static Function<String, ChrBaseRegion> factory()
    {
        return s -> new BedLine(s).toRegion();
    }

    private final String[] data;

    public BedLine(String line)
    {
        data = line.split(TSV_DELIM, -1);

        if(data.length < 3)
        {
            throw new IllegalArgumentException("invalid slice BED entry: " + line);
        }
    }

    public Chromosome chromosome()
    {
        if(!HumanChromosome.contains(data[0]))
        {
            return null;
        }

        return HumanChromosome.fromString(data[0]);
    }

    public int start()
    {
        return Integer.parseInt(data[1]) + 1; // as per convention
    }

    public int end()
    {
        return Integer.parseInt(data[2]);
    }

    public ChrBaseRegion toRegion()
    {
        Chromosome chromosome = chromosome();
        if(chromosome == null)
        {
            return null;
        }
        return new ChrBaseRegion(chromosome.toString(), start(), end());
    }

    public String dataValueOrBlank(int position)
    {
        if(position < 0 || position >= data.length)
        {
            return "";
        }
        return data[position];
    }
}
