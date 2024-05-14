package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

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

    public static SglMapping from(final String mappingStr, byte breakendOrientation)
    {
        final String[] items = mappingStr.split(INS_SEQ_DATA_DELIM, -1);

        if(items.length != 4)
            return null;

        final String[] location = items[0].split(":");

        if(location.length != 2)
            return null;

        if(items[1].isEmpty())
            return null;

        String chromosome = location[0];

        if(!HumanChromosome.contains(chromosome))
            return null;

        int position = Integer.parseInt(location[1]);

        final String cigar = items[2];

        byte strandSequence = items[1].equals("+") ? POS_ORIENT : NEG_ORIENT;
        byte orientation = (strandSequence == breakendOrientation) ? NEG_ORIENT : POS_ORIENT;

        if(orientation == POS_ORIENT)
            position += calcCigarAlignedLength(cigar);

        int qualScore = !items[3].isEmpty() ? Integer.parseInt(items[3]) : 0;

        return new SglMapping(chromosome, position, orientation, cigar, qualScore);
    }

    public static void convertFromInsertSequenceAlignments(
            final List<SglMapping> mappings, final String alignments, byte breakendOrientation)
    {
        if(alignments.isEmpty())
            return;

        final String[] mappingStrings = alignments.split(INS_SEQ_MAPPING_DELIM, -1);

        for(final String mappingData : mappingStrings)
        {
            final SglMapping mapping = from(mappingData, breakendOrientation);

            if(mapping != null)
                mappings.add(mapping);
        }
    }

    public boolean possibleLink(final SglMapping other)
    {
        return isPossibleLink(
                Chromosome, Position, Orientation, other.Chromosome, other.Position, other.Orientation, MIN_TEMPLATED_INSERTION_LENGTH);
    }

    public String toString()
    {
        return String.format("%s:%d:%d qual=%d", Chromosome, Position, Orientation, QualScore);
    }

}
