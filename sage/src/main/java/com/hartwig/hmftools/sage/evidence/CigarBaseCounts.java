package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

class CigarBaseCounts
{
    public int SoftClipStart;
    public int SoftClipEnd;
    public int AdjustedBases;
    public int AlignedBases;

    public CigarBaseCounts(final Cigar cigar)
    {
        SoftClipStart = 0;
        SoftClipEnd = 0;
        AdjustedBases = 0;
        AlignedBases = 0;

        for(int i = 0; i < cigar.getCigarElements().size(); ++i)
        {
            CigarElement element = cigar.getCigarElements().get(i);

            switch(element.getOperator())
            {
                case S:
                {
                    if(i == 0)
                        SoftClipStart += element.getLength();
                    else
                        SoftClipEnd += element.getLength();
                    break;
                }

                case D:
                case N:
                    AdjustedBases -= element.getLength();
                    AlignedBases += element.getLength();
                    break;

                case M:
                    AlignedBases += element.getLength();
                    break;

                case I:
                    AdjustedBases += element.getLength();

                default:
                    break;
            }
        }
    }

    public String toString()
    {
        return format("sc(%d - %d) adjusted(%d) aligned(%d)", SoftClipStart, SoftClipEnd, AdjustedBases, AlignedBases);
    }
}
