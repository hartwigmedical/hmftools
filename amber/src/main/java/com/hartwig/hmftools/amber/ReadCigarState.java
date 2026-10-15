package com.hartwig.hmftools.amber;

import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.CigarElement;

public class ReadCigarState
{
    public final List<CigarElement> Elements;
    public final int BaseLength;
    public int ReadIndex;
    public int RefPosition;
    public int CigarIndex; // cigar element
    public int ElementIndex; // index within cigar element

    public ReadCigarState(final List<CigarElement> cigarElements, final int readStart, final int baseLength)
    {
        Elements = cigarElements;
        BaseLength = baseLength;
        RefPosition = readStart;
        ReadIndex = 0;
        CigarIndex = 0;
        ElementIndex = 0;

        if(cigarElements.get(0).getOperator() == S)
        {
            CigarIndex = 1;
            ReadIndex = cigarElements.get(0).getLength();
        }
    }

    public void moveToPosition(int position)
    {
        if(position <= RefPosition)
            return;

        for(; CigarIndex < Elements.size(); ++CigarIndex)
        {
            CigarElement element = Elements.get(CigarIndex);

            // check case for end of an aligned section
            if(ElementIndex == element.getLength() - 1)
            {
                ElementIndex = 0;

                if(element.getOperator().consumesReadBases())
                    ++ReadIndex;

                if(element.getOperator().consumesReferenceBases())
                    ++RefPosition;

                // return if the required position falls into a deleted section
                if(CigarIndex < Elements.size() - 1)
                {
                    CigarElement nextElement = Elements.get(CigarIndex + 1);

                    if(nextElement.getOperator() == D)
                    {
                        if(RefPosition == position)
                        {
                            ++CigarIndex;
                            return;
                        }
                    }
                }

                continue;
            }

            if(!element.getOperator().consumesReferenceBases())
            {
                ReadIndex += element.getLength();
                ElementIndex = 0;
                continue;
            }

            for(; ElementIndex < element.getLength(); ++ElementIndex)
            {
                if(RefPosition == position)
                    break;

                ++RefPosition;

                if(element.getOperator().consumesReadBases())
                    ++ReadIndex;
            }

            if(ReadIndex >= BaseLength)
                return;

            if(RefPosition == position)
                return;

            ElementIndex = 0;
        }
    }

    public boolean isValid() { return CigarIndex >= 0 && CigarIndex < Elements.size() && ReadIndex < BaseLength; }

    public boolean isDelete()
    {
        if(!isValid())
            return false;

        return Elements.get(CigarIndex).getOperator() == D;
    }

    public boolean isInsert()
    {
        if(!isValid())
            return false;

        // must be at the end of an aligned section with an insert following
        CigarElement element = Elements.get(CigarIndex);
        if(element.getOperator() == M && ElementIndex == element.getLength() - 1 && CigarIndex < Elements.size() - 1)
            return Elements.get(CigarIndex + 1).getOperator() == I;

        return false;
    }

    public boolean isIndel() { return isDelete() || isInsert(); }

    public String toString()
    {
        CigarElement element = isValid() ? Elements.get(CigarIndex) : null;

        return format("cigar(%d of %d) element(%d of %s) readIndex(%d) refPos(%d)",
                CigarIndex, Elements.size(), ElementIndex, element, ReadIndex, RefPosition);
    }
}
