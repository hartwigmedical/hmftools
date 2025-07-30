package com.hartwig.hmftools.common.region;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.Function;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class TaggedRegion extends BaseRegion
{
    private static class TaggedChromosomeRegion extends ChrBaseRegion
    {
        String tag;

        public TaggedChromosomeRegion(final Chromosome chromosome, final int posStart, final int posEnd, final String tag)
        {
            super(chromosome.toString(), posStart, posEnd);
            this.tag = tag;
        }

        TaggedRegion toTaggedRegion()
        {
            return new TaggedRegion(start(), end(), tag);
        }
    }

    public static Map<Chromosome, List<TaggedRegion>> loadRegionsFromBedFile(String bedFile)
    {
        Function<String, TaggedChromosomeRegion> factory = s ->
        {
            BedLine line = new BedLine(s);
            Chromosome chromosome = line.chromosome();
            if(chromosome == null)
            {
                return null;
            }
            return new TaggedChromosomeRegion(chromosome, line.start(), line.end(), line.dataValueOrBlank(3));
        };
        Function<TaggedChromosomeRegion, TaggedRegion> converter = TaggedChromosomeRegion::toTaggedRegion;
        return BedFileReader.loadBedFileChrMap(bedFile, true, factory, converter);
    }

    public final String mTag;

    public TaggedRegion(final int posStart, final int posEnd, final String mTag)
    {
        super(posStart, posEnd);
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
