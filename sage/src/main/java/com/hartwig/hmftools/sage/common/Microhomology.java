package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import htsjdk.samtools.SAMRecord;

public class Microhomology
{
    public final String Bases;
    public final int Length;

    public Microhomology(final String bases, final int length)
    {
        Bases = bases;
        Length = length;
    }

    public String toString()
    {
        return format("%s length(%d)", Bases, Length);
    }

    public static Microhomology findHomology(final SimpleVariant variant, final SAMRecord read, int varReadIndex)
    {
        if(!variant.isIndel())
            return null;

        int indelAltLength = abs(variant.indelLength());

        StringBuilder homology = null;
        String indelBases = variant.isInsert() ? variant.alt().substring(1) : variant.ref().substring(1);

        final byte[] readBases = read.getReadBases();

        // start looking in the read in the first base after the variant
        int homReadIndexStart = variant.isInsert() ? varReadIndex + indelAltLength + 1 : varReadIndex + 1;
        int homReadIndex = homReadIndexStart;

        for(int i = 0; i < indelBases.length() && homReadIndex < readBases.length; ++i, ++homReadIndex)
        {
            byte indelBase = (byte)indelBases.charAt(i);
            byte postIndelReadBase = readBases[homReadIndex];

            if(indelBase != postIndelReadBase)
                break;

            if(homology == null)
                homology = new StringBuilder();

            homology.append((char)indelBase);
        }

        if(homology == null)
            return null;

        String homologyBases = homology.toString();
        int homologyLength = homologyBases.length();

        if(homologyLength == indelAltLength && homReadIndex < readBases.length)
        {
            // continue searching for repeats of the homology to find the point of right-alignment, allow partials at the end
            boolean matched = true;

            while(matched)
            {
                int matchCount = 0;

                for(int i = 0; i < indelAltLength && homReadIndex < readBases.length; ++i, ++homReadIndex)
                {
                    byte homologyBase = (byte)homologyBases.charAt(i);
                    byte postIndelReadBase = readBases[homReadIndex];

                    if(homologyBase != postIndelReadBase)
                    {
                        matched = false;
                        break;
                    }

                    ++matchCount;
                }

                if(matchCount > 0)
                    homologyLength += matchCount;
                else
                    break;
            }
        }

        return new Microhomology(homologyBases, homologyLength);
    }

}
