package com.hartwig.hmftools.bamtools.remapper;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class AltAlignment
{
    // consider moving into common and using in Linx instead of its SglMapping class
    public final String Chromosome;
    public final int Position;
    public final Orientation Orient;
    public final String Cigar;
    public final int EditDistance;

    public AltAlignment(final String chromosome, final int position, final Orientation orientation, final String cigar, final int editDistance)
    {
        Chromosome = chromosome;
        Position = position;
        Orient = orientation;
        Cigar = cigar;
        EditDistance = editDistance;
    }

    private static final String ITEM_DELIM = ",";

    public static List<AltAlignment> fromLocationTag(final String locationTag)
    {
        if(locationTag.isEmpty())
            return Collections.emptyList();

        final String[] mappingStrings = locationTag.split(";", -1);

        List<AltAlignment> alignments = Lists.newArrayList();

        for(final String mappingData : mappingStrings)
        {
            AltAlignment mapping = fromBwaString(mappingData);

            if(mapping != null)
                alignments.add(mapping);
        }

        return alignments;
    }

    private static AltAlignment fromBwaString(final String mappingStr)
    {
        // expected BWA format: // 16,+24008715,44S28M46S,0;X,+133232624,44S27M47S,0;12,+54042138,37S35M46S,2;4,-84437081,46S25M47S,0;
        final String[] items = mappingStr.split(ITEM_DELIM, 4);

        if(items.length != 4)
            return null;

        String chromosome = items[0];

        if(!HumanChromosome.contains(chromosome))
            return null;

        String orientPos = items[1];
        Orientation orientation = Orientation.fromChar(orientPos.charAt(0));

        int position = Integer.parseInt(orientPos.substring(1));

        final String cigar = items[2];

        int mapQual = Integer.parseInt(items[3]);

        return new AltAlignment(chromosome, position, orientation, cigar, mapQual);
    }

    public String toString()
    {
        return format("%s:%d:%d mq=%d", Chromosome, Position, Orient.asByte(), EditDistance);
    }
}
