package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

public class SglMapping
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String Cigar;
    public final int QualScore;

    public SglMapping(final String chromosome, final int position, final byte orientation, final String cigar, final int qualScore)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        Cigar = cigar;
        QualScore = qualScore;
    }

    private static final String INS_SEQ_MAPPING_DELIM = ",";
    private static final String INS_SEQ_DATA_DELIM = "\\|";

    // expected format: 4:9973661|-|26S37M11S|19,1:184636186|+|37M37S|

    private static SglMapping from(final String mappingStr)
    {
        final String[] items = mappingStr.split(INS_SEQ_DATA_DELIM, -1);

        if(items.length != 4)
            return null;

        final String[] location = items[0].split(":");

        if(location.length != 2)
            return null;

        if(items[1].isEmpty() || items[3].isEmpty())
            return null;

        String chromosome = location[0];
        int position = Integer.parseInt(location[1]);
        byte orientation = items[1].equals("+") ? POS_ORIENT : NEG_ORIENT;
        String cigar = items[2];
        int qualScore = Integer.parseInt(items[3]);

        return new SglMapping(chromosome, position, orientation, cigar, qualScore);
    }

    public static void convertFromInsertSequenceAlignments(final List<SglMapping> mappings, final String alignments)
    {
        if(alignments.isEmpty())
            return;

        final String[] mappingStrings = alignments.split(INS_SEQ_MAPPING_DELIM, -1);

        for(final String mappingData : mappingStrings)
        {
            final SglMapping mapping = from(mappingData);

            if(mapping != null)
                mappings.add(mapping);
        }
    }

    public boolean possibleLink(final SglMapping other)
    {
        return isPossibleLink(
                Chromosome, Position, Orientation, other.Chromosome, other.Position, other.Orientation, MIN_TEMPLATED_INSERTION_LENGTH);
    }
}
