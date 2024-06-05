package com.hartwig.hmftools.gripss.rm;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.calcCigarAlignedLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipped;

import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;

public class AlignmentData
{
    public final String Chromosome;
    public final BaseRegion Region;
    public final String CigarStr;

    private static final String INS_SEQ_MAPPING_DELIM = ",";
    private static final String INS_SEQ_MAPPING_DELIM_2 = ";";
    private static final String INS_SEQ_DATA_DELIM = "\\|";

    public AlignmentData(final String chromosome, final BaseRegion region, final String cigar)
    {
        Chromosome = chromosome;
        Region = region;
        CigarStr = cigar;
    }

    public int overlappingBases(final BaseRegion region)
    {
        return min(Region.end(), region.end()) - max(Region.start(), region.start()) + 1;
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
        int cigarLength = calcCigarAlignedLength(cigar);

        return new AlignmentData(chromosome, new BaseRegion(position, position + cigarLength - 1), cigar);
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

    public String extractMatchedBases(final String insertSequence)
    {
        Cigar cigar = CigarUtils.cigarFromStr(CigarStr);

        int matchStartPos = leftSoftClipped(cigar) ? cigar.getFirstCigarElement().getLength() : 0;
        int matchBases = cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() == M || x.getOperator() == I)
                .mapToInt(x -> x.getLength()).sum();

        int insSeqLength = insertSequence.length();
        int startPos = min(matchStartPos, insSeqLength);
        int endPos = min(matchStartPos + matchBases, insSeqLength);

        return insertSequence.substring(startPos, endPos);
    }
}
