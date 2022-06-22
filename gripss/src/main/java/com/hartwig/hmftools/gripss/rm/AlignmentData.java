package com.hartwig.hmftools.gripss.rm;

import static com.hartwig.hmftools.common.samtools.CigarTraversal.calcCigarLength;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class AlignmentData
{
    public final String Chromosome;
    public final BaseRegion Region;

    private static final String INS_SEQ_MAPPING_DELIM = ",";
    private static final String INS_SEQ_MAPPING_DELIM_2 = ";";
    private static final String INS_SEQ_DATA_DELIM = "\\|";

    public AlignmentData(final String chromosome, final BaseRegion region)
    {
        Chromosome = chromosome;
        Region = region;
    }

    public static AlignmentData fromAlignment(final String alignment)
    {
        // expected format: 4:9973661|-|26S37M11S|19,1:184636186|+|37M37S|
        final String[] values = alignment.split(INS_SEQ_DATA_DELIM, -1);

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

        final String cigar = values[2];
        int cigarLength = calcCigarLength(cigar);

        return new AlignmentData(chromosome, new BaseRegion(position, position + cigarLength - 1));
    }

    public static List<AlignmentData> fromInsertSequenceAlignments(final String alignmentsStr)
    {
        final List<AlignmentData> alignments = Lists.newArrayList();

        if(alignmentsStr == null || alignmentsStr.isEmpty())
            return alignments;

        final String[] mappingStrings = alignmentsStr.contains(INS_SEQ_MAPPING_DELIM) ?
                alignmentsStr.split(INS_SEQ_MAPPING_DELIM, -1) : alignmentsStr.split(INS_SEQ_MAPPING_DELIM_2, -1);

        for(final String mappingData : mappingStrings)
        {
            final AlignmentData alignment = fromAlignment(mappingData);

            if(alignment != null)
                alignments.add(alignment);
        }

        return alignments;
    }

    public String toString()
    {
        return String.format("%s:%s", Chromosome, Region.toString());
    }

}
