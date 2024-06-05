package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.Cigar;

public class AlternativeAlignment
{
    // consider moving into common and using in Linx instead of its SglMapping class
    public final String Chromosome;
    public final int Position;
    public final Orientation Orient;
    public final String Cigar;
    public final int MapQual;

    private Cigar mCigar;

    public AlternativeAlignment(final String chromosome, final int position, final Orientation orientation, final String cigar, final int mapQual)
    {
        Chromosome = chromosome;
        Position = position;
        Orient = orientation;
        Cigar = cigar;
        MapQual = mapQual;

        mCigar = null;
    }

    private static final String MAPPING_DELIM = ";";
    private static final String ITEM_DELIM = ",";
    private static final String SUB_ITEM_DELIM = "|";

    public Cigar cigar()
    {
        if(mCigar == null)
            mCigar = CigarUtils.cigarFromStr(Cigar);

        return mCigar;
    }

    private static AlternativeAlignment fromBwaString(final String mappingStr)
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

        return new AlternativeAlignment(chromosome, position, orientation, cigar, mapQual);
    }

    public static List<AlternativeAlignment> fromLocationTag(final String locationTag)
    {
        if(locationTag.isEmpty())
            return Collections.emptyList();

        final String[] mappingStrings = locationTag.split(MAPPING_DELIM, -1);

        List<AlternativeAlignment> alignments = Lists.newArrayList();

        for(final String mappingData : mappingStrings)
        {
            AlternativeAlignment mapping = fromBwaString(mappingData);

            if(mapping != null)
                alignments.add(mapping);
        }

        return alignments;
    }

    public String toString()
    {
        return format("%s:%d:%d mq=%d", Chromosome, Position, Orient, MapQual);
    }

    public String vcfString()
    {
        // in form expected by Linx and other downstream components: 4:9973661|-|26S37M11S|19
        StringJoiner sj = new StringJoiner(SUB_ITEM_DELIM);
        sj.add(format("%s:%d", Chromosome, Position));
        sj.add(String.valueOf(Orient.asChar()));
        sj.add(Cigar);
        sj.add(String.valueOf(MapQual));
        return sj.toString();
    }

    public static String toVcfTag(final List<AlternativeAlignment> alts)
    {
        if(alts == null || alts.isEmpty())
            return "";

        return alts.stream().map(x -> x.vcfString()).collect(Collectors.joining(ITEM_DELIM));
    }

    public static List<AlternativeAlignment> fromVcfTag(final String vcfTag)
    {
        if(vcfTag.isEmpty())
            return Collections.emptyList();

        final String[] altAlignmentStr = vcfTag.split(ITEM_DELIM, -1);

        List<AlternativeAlignment> alignments = Lists.newArrayList();

        for(final String altAlignStr : altAlignmentStr)
        {
            AlternativeAlignment altAlignment = fromVcfString(altAlignStr);

            if(altAlignment != null)
                alignments.add(altAlignment);
        }

        return alignments;
    }

    private static AlternativeAlignment fromVcfString(final String alignmentStr)
    {
        // expected format: 4:9973661|-|26S37M11S|19,1:184636186|+|37M37S|
        final String[] values = alignmentStr.split("\\|", -1);

        if(values.length != 4)
            return null;

        final String[] location = values[0].split(":");

        if(location.length != 2)
            return null;

        if(values[1].isEmpty())
            return null;

        String chromosome = location[0];

        if(!HumanChromosome.contains(chromosome))
            return null;

        int position = Integer.parseInt(location[1]);

        Orientation orientation = Orientation.fromChar(values[1].charAt(0));

        final String cigar = values[2];

        int mapQual = Integer.parseInt(values[3]);

        return new AlternativeAlignment(chromosome, position, orientation, cigar, mapQual);
    }

    public boolean matches(final AlternativeAlignment other)
    {
        return Chromosome.equals(other.Chromosome) && Position == other.Position && Orient == other.Orient
                && Cigar.equals(other.Cigar) && MapQual == other.MapQual;
    }
}
