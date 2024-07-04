package com.hartwig.hmftools.sage.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

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
        return format("%s length=%d", Bases, Length);
    }

    public static int findLeftHomologyShift(
            final SimpleVariant variant, final RefSequence refSequence, final byte[] readBases, int varReadIndex)
    {
        if(!variant.isIndel())
            return 0;

        int leftAlignOffset = 0;

        int currentPosition = variant.Position;

        byte lastRefBase = (byte)variant.Ref.charAt(variant.refLength() - 1);

        String altBases = variant.Alt;
        byte lastAltBase = (byte)altBases.charAt(variant.altLength() - 1);

        while(lastRefBase == lastAltBase)
        {
            --currentPosition;
            ++leftAlignOffset;

            if(currentPosition < refSequence.Start)
                break;

            lastRefBase = refSequence.base(currentPosition);

            if(varReadIndex - leftAlignOffset < 0)
                break;

            char prevAltBase = (char)readBases[varReadIndex - leftAlignOffset];
            altBases = prevAltBase + altBases.substring(0, altBases.length() - 1);
            lastAltBase = (byte)altBases.charAt(variant.altLength() - 1);
        }

        return leftAlignOffset;
    }

    public static Microhomology findHomology(final SimpleVariant variant, final byte[] readBases, int varReadIndex, boolean applyIndelLength)
    {
        if(!variant.isIndel())
            return null;

        int indelAltLength = variant.indelLengthAbs();

        StringBuilder homology = null;
        String indelBases = variant.isInsert() ? variant.alt().substring(1) : variant.ref().substring(1);

        // start looking in the read in the first base after the variant
        int homReadIndexStart = varReadIndex + 1;

        if(applyIndelLength && variant.isInsert())
            homReadIndexStart += indelAltLength;

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
