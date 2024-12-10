package com.hartwig.hmftools.common.sequencing;

public class BiomodalBamUtils
{
    public static int LOW_QUAL_CUTOFF = 30;

    public static String encodeMMTag(final StringBuilder readStr, byte modCBase)
    {
        StringBuilder MMTag = new StringBuilder("MM:Z:C+C.");
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
        return MMTag.toString();
    }

    public static void decodeMMTag(final StringBuilder readStr, final String mmTag, byte modCBase)
    {

    }
}
