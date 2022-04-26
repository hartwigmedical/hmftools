package com.hartwig.hmftools.linx.types;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.chromosomeRank;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

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
        byte orientation = items[1].equals("+") ? POS_ORIENT : NEG_ORIENT;
        final String cigar = items[2];

        if(orientation != breakendOrientation)
            position += calcCigarLength(cigar);

        int qualScore = !items[3].isEmpty() ? Integer.parseInt(items[3]) : 0;

        return new SglMapping(chromosome, position, orientation, cigar, qualScore);
    }

    public static int calcCigarLength(final String cigarStr)
    {
        int index = 0;
        int baseLength = 0;
        String basesStr = "";
        while(index < cigarStr.length())
        {
            char c = cigarStr.charAt(index);
            boolean isAddItem = (c == 'D' || c == 'M');
            boolean isIgnoreItem = (c == 'I' || c == 'N' || c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X');

            if(isAddItem)
            {
                try { baseLength += Integer.parseInt(basesStr); } catch (Exception e) {}
                basesStr = "";
            }
            else if(isIgnoreItem)
            {
                basesStr = "";
            }
            else
            {
                basesStr += c;
            }

            ++index;
        }

        return baseLength;
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
