package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.D;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;

import htsjdk.samtools.CigarElement;

public class IndelCoords
{
    public final int PosStart;
    public final int PosEnd;
    public final int Length;

    private String mInsertedBases;

    public IndelCoords(final int posStart, final int posEnd, final int length)
    {
        PosStart = posStart;
        PosEnd = posEnd;
        Length = length;
        mInsertedBases = null;
    }

    public boolean isInsert() { return PosEnd == PosStart + 1; }
    public boolean isDelete() { return !isInsert(); }

    public String insertedBases() { return mInsertedBases != null ? mInsertedBases : ""; }
    public void setInsertedBases(final String bases) { mInsertedBases = bases; }

    public boolean matchesJunction(int position, Orientation orientation)
    {
        return orientation.isForward() ? PosStart == position : PosEnd == position;
    }

    public boolean matches(final IndelCoords other)
    {
        return PosStart == other.PosStart && PosEnd == other.PosEnd && Length == other.Length;
    }

    public String toString()
    {
        return format("%s(%d - %d) len(%d)", isDelete() ? "delete" : "insert", PosStart, PosEnd, Length);
    }


    public static IndelCoords findIndelCoords(final int readStart, final List<CigarElement> cigarElements, int minIndelLength)
    {
        int maxIndelLength = 0;

        // find the location of the internal delete or insert matching the max indel length
        int indelStartPos = readStart - 1;
        int indelEndPos = 0;

        int refPosition = readStart;

        for(CigarElement element : cigarElements)
        {
            if(element.getOperator().isIndel())
            {
                if(element.getLength() > maxIndelLength)
                {
                    maxIndelLength = element.getLength();

                    // indel start is the last ref base, indel end is the next ref base
                    indelStartPos = refPosition - 1;

                    if(element.getOperator() == D)
                        indelEndPos = indelStartPos + element.getLength() + 1;
                    else
                        indelEndPos = indelStartPos + 1;
                }
            }

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();
        }

        if(indelEndPos <= indelStartPos || maxIndelLength < minIndelLength)
            return null;

        return new IndelCoords(indelStartPos, indelEndPos, maxIndelLength);
    }

    public static IndelCoords findIndelCoords(final int readStart, final List<CigarElement> cigarElements, final CigarElement specificIndel)
    {
        // find the location of the internal delete or insert matching the max indel length
        int indelStartPos = readStart - 1;
        int indelEndPos = 0;
        int matchedIndelLength = 0;

        int refPosition = readStart;

        for(CigarElement element : cigarElements)
        {
            if(element.getOperator() == specificIndel.getOperator() && abs(specificIndel.getLength() - element.getLength()) <= 1)
            {
                // indel start is the last ref base, indel end is the next ref base
                indelStartPos = refPosition - 1;

                if(element.getOperator() == D)
                    indelEndPos = indelStartPos + element.getLength() + 1;
                else
                    indelEndPos = indelStartPos + 1;

                matchedIndelLength = element.getLength();
                break;
            }

            if(element.getOperator().consumesReferenceBases())
                refPosition += element.getLength();
        }

        return new IndelCoords(indelStartPos, indelEndPos, matchedIndelLength);
    }
}
