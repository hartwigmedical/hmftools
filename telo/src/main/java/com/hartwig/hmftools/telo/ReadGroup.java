package com.hartwig.hmftools.telo;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    public final List<SAMRecord> Reads = new ArrayList<>();

    public ReadGroup(final SAMRecord read)
    {
        Reads.add(read);
    }

    public final String id() { return Reads.get(0).getReadName(); }

    public boolean isComplete()
    {
        return !Reads.get(0).getReadPairedFlag() || Reads.size() == 3 || (Reads.size() == 2 && !hasSuppAlignment(Reads));
    }

    public static boolean hasSuppAlignment(final List<SAMRecord> reads)
    {
        return reads.stream().anyMatch(x -> x.getSupplementaryAlignmentFlag());
    }

    //public boolean isDuplicate() { return Reads.stream().anyMatch(x -> x.isDuplicate()); }

    public String toString() { return String.format("%s reads(%d) complete(%s)", id(), Reads.size(), isComplete()); }

    /*
    public String findOtherChromosome(final String chromosome)
    {
        for(ReadRecord read : Reads)
        {
            if(!read.mateChromosome().equals(chromosome))
                return read.mateChromosome();

            if(read.hasSuppAlignment())
                return suppAlignmentChromosome(read.getSuppAlignment());
        }

        return null;
    }*/

    public static final String SUPP_ALIGNMENT_DELIM = ",";

    public static String suppAlignmentChromosome(final String suppAlignment)
    {
        if(suppAlignment == null)
            return null;

        final String[] items = suppAlignment.split(SUPP_ALIGNMENT_DELIM);
        return items.length >= 5 ? items[0] : null;
    }
}

