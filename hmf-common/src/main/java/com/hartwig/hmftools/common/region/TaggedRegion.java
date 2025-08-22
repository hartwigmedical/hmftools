package com.hartwig.hmftools.common.region;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;

import com.hartwig.hmftools.common.genome.chromosome.Chromosomal;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TaggedRegion extends ChrBaseRegion implements Chromosomal
{
    private static final Logger LOGGER = LogManager.getLogger(TaggedRegion.class);

    public static Map<Chromosome, List<TaggedRegion>> loadRegionsFromBedFile(final String bedFile)
    {
        Function<String, TaggedRegion> factory = s ->
        {
            BedLine line = new BedLine(s);
            Chromosome chromosome = line.chromosome();
            if(chromosome == null)
            {
                return null;
            }
            return new TaggedRegion(chromosome.toString(), line.start(), line.end(), line.dataValueOrBlank(3));
        };
        try
        {
            List<TaggedRegion> fromBed = BedFileReader.loadBedFile(bedFile, true, factory);
            return Chromosomal.toPerChromosomeLists(fromBed);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load BED file({}): {}", bedFile, e.toString());
            return null;
        }
    }

    public final String mTag;

    public TaggedRegion(String chromosome, final int posStart, final int posEnd, final String mTag)
    {
        super(chromosome, posStart, posEnd);
        this.mTag = mTag;
    }

    public String formatted()
    {
        return String.format("%s:%d-%d", mTag, start(), end());
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        if(!super.equals(o))
        {
            return false;
        }
        final TaggedRegion that = (TaggedRegion) o;
        return Objects.equals(mTag, that.mTag);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(super.hashCode(), mTag);
    }
}
