package com.hartwig.hmftools.lilac.misc;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import htsjdk.samtools.SAMRecord;

public class LilacTestUtils
{
    public static SAMRecord buildSamRecord(int alignmentStart, String cigar, String readString, String qualities)
    {
        SAMRecord record = new SAMRecord(null);
        record.setReadName("READNAME");
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        return record;
    }

    public static HlaSequenceLoci createSequenceLoci(final HlaSequence sequences)
    {
        return HlaSequenceLoci.create(sequences.Allele, sequences.getRawSequence(), sequences.getRawSequence());
    }
    public static String buildTargetSequence(final String sequence, final List<Integer> indices)
    {
        if(indices.size() >= sequence.length())
            return "";

        StringBuilder targetSeq = new StringBuilder();
        indices.forEach(x -> targetSeq.append(sequence.charAt(x)));
        return targetSeq.toString();
    }

    public static List<Integer> buildLoci(final String sequence)
    {
        List<Integer> indices = Lists.newArrayList();

        for(int i = 0; i < sequence.length(); ++i)
        {
            indices.add(i);
        }

        return indices;
    }


    public static List<String> buildTargetSequences(final String sequence, final List<Integer> indices)
    {
        List<String> sequences = Lists.newArrayList();
        if(indices.size() > sequence.length())
            return sequences;

        indices.forEach(x -> sequences.add(String.valueOf(sequence.charAt(x))));
        return sequences;
    }

    public static List<Integer> buildTargetLoci(final String sequence, final String refSequence)
    {
        List<Integer> indices = Lists.newArrayList();

        if(sequence.length() != refSequence.length())
            return indices;

        for(int i = 0; i < sequence.length(); ++i)
        {
            if(sequence.charAt(i) != refSequence.charAt(i))
            {
                indices.add(i);
            }
        }

        return indices;
    }
}
