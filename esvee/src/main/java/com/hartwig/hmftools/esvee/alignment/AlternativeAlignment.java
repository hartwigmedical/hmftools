package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class AlternativeAlignment
{
    // consider moving into common and using in Linx instead of its SglMapping class
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String Cigar;
    public final int MapQual;

        public AlternativeAlignment(final String chromosome, final int position, final byte orientation, final String cigar, final int mapQual)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        Cigar = cigar;
        MapQual = mapQual;
    }

    private static final String MAPPING_DELIM = ";";
    private static final String ITEM_DELIM = ",";
    private static final char POS_ORIENT_CHAR = '+';
    private static final char NEG_ORIENT_CHAR = '-';

    // expected BWA format: // 16,+24008715,44S28M46S,0;X,+133232624,44S27M47S,0;12,+54042138,37S35M46S,2;4,-84437081,46S25M47S,0;

    public static AlternativeAlignment from(final String mappingStr)
    {
        final String[] items = mappingStr.split(ITEM_DELIM, 4);

        if(items.length != 4)
            return null;

        String chromosome = items[0];

        if(!HumanChromosome.contains(chromosome))
            return null;

        String orientPos = items[1];
        byte orientation = orientPos.charAt(0) == POS_ORIENT_CHAR ? POS_ORIENT : NEG_ORIENT;

        int position = Integer.parseInt(orientPos.substring(1));

        final String cigar = items[2];

        int mapQual = Integer.parseInt(items[3]);

        return new AlternativeAlignment(chromosome, position, orientation, cigar, mapQual);
    }

    public static List<AlternativeAlignment> convertLocationTag(final String locationTag)
    {
        if(locationTag.isEmpty())
            return Collections.emptyList();

        final String[] mappingStrings = locationTag.split(MAPPING_DELIM, -1);

        List<AlternativeAlignment> alignments = Lists.newArrayList();

        for(final String mappingData : mappingStrings)
        {
            AlternativeAlignment mapping = from(mappingData);

            if(mapping != null)
                alignments.add(mapping);
        }

        return alignments;
    }

    public String toString()
    {
        return format("%s:%d:%d mq=%d", Chromosome, Position, Orientation, MapQual);
    }

    public String vcfString()
    {
        // in form expected by Linx and other downstream components: 4:9973661|-|26S37M11S|19
        // make common class for this
        StringJoiner sj = new StringJoiner("|");
        sj.add(format("%s:%d", Chromosome, Position));
        sj.add(String.valueOf(Orientation == POS_ORIENT ? POS_ORIENT_CHAR : NEG_ORIENT_CHAR));
        sj.add(Cigar);
        sj.add(String.valueOf(MapQual));
        return sj.toString();
    }

    public static String altAlignmentsStr(final List<AlternativeAlignment> alts)
    {
        if(alts == null || alts.isEmpty())
            return "";

        return alts.stream().map(x -> x.vcfString()).collect(Collectors.joining(ITEM_DELIM));
    }
}
