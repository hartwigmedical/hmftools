package com.hartwig.hmftools.common.sequencing;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import org.apache.commons.lang3.NotImplementedException;

import htsjdk.samtools.SAMRecord;

public class BiomodalBamUtils
{
    public static int LOW_QUAL_CUTOFF = 30;
    public static byte MODC_BASE = (byte) 'X';

    private static String MM_TAG_PREFIX = "MM:Z:C+C.";

    public static String encodeMMTag(final StringBuilder readStr, boolean isForward, byte modCBase)
    {
        if(!isForward)
            readStr.reverse();

        StringBuilder MMTag = new StringBuilder(MM_TAG_PREFIX);
        int skip = 0;
        for(int i = 0; i < readStr.length(); i++)
        {
            char base = readStr.charAt(i);
            if(base == 'C')
            {
                skip++;
            }
            else if((byte) base == modCBase)
            {
                readStr.setCharAt(i, 'C');
                MMTag.append(',');
                MMTag.append(skip);
                skip = 0;
            }
        }

        MMTag.append(';');

        if(!isForward)
            readStr.reverse();

        return MMTag.toString();
    }

    public static byte[] decodeMMTag(final SAMRecord read, byte modCBase)
    {
        // TODO:
        if(true)
        {
            throw new NotImplementedException("TODO");
        }

        String mmTag = read.getStringAttribute("MM");
        if(mmTag == null)
        {
            return read.getReadBases();
        }

        if(mmTag.length() == MM_TAG_PREFIX.length() + 1)
        {
            return read.getReadBases();
        }

        String[] mmSkipStrs = mmTag.substring(MM_TAG_PREFIX.length() + 1, mmTag.length() - 1).split("[,]");
        int[] mmSkips = new int[mmSkipStrs.length];
        for(int i = 0; i < mmSkipStrs.length; i++)
            mmSkips[i] = Integer.parseInt(mmSkipStrs[i]);

        boolean isForward = !read.getReadNegativeStrandFlag();
        byte[] readBases = read.getReadBases();
        byte[] decodedBases = new byte[readBases.length];

        int i = isForward ? 0 : readBases.length - 1;
        int inc = isForward ? 1 : -1;
        int skipIdx = 0;
        for(; i >= 0 && i < readBases.length; i += inc)
        {
            if(skipIdx >= mmSkips.length)
            {
                decodedBases[i] = readBases[i];
                continue;
            }

            byte base = readBases[i];
            if(!isForward)
                base = swapDnaBase(base);

            if(base != (byte) 'C')
            {
                decodedBases[i] = readBases[i];
                continue;
            }

            if(mmSkips[skipIdx] > 0)
            {
                mmSkips[skipIdx]--;
                continue;
            }
        }

        return decodedBases;
    }
}
