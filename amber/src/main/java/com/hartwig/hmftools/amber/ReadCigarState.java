package com.hartwig.hmftools.amber;

import static java.lang.Math.min;
import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.CigarElement;

public class ReadCigarState
{
    public final List<CigarElement> Elements; // left null if a single aligned section
    public final int BaseLength;
    public int ReadIndex;
    public int RefPosition;
    public int CigarIndex; // cigar element
    public int ElementIndex; // index within cigar element

    public ReadCigarState(final List<CigarElement> cigarElements, final int readStart, final int baseLength)
    {
        if(cigarElements.size() == 1 && cigarElements.get(0).getOperator() == M)
        {
            Elements = null;
        }
        else
        {
            Elements = cigarElements;
        }

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

        if(Elements == null)
        {
            moveToPositionAligned(position);
            return;
        }

        for(; CigarIndex < Elements.size(); ++CigarIndex)
        {
            CigarElement element = Elements.get(CigarIndex);

            if(!element.getOperator().consumesReferenceBases())
            {
                ReadIndex += element.getLength();
                ElementIndex = -1;
                continue;
            }

            // eg current pos = 10, required pos = 15, posDiff = 5
            // current element index = 2, length = 7, remaining length = 4
            int posDifference = position - RefPosition;
            int remainingElement = element.getLength() - ElementIndex - 1;
            int posMove = min(posDifference, remainingElement);

            RefPosition += posMove;

            if(element.getOperator().consumesReadBases())
                ReadIndex += posMove;

            ElementIndex += posMove;

            if(ReadIndex >= BaseLength)
                return;

            if(RefPosition == position)
                return;

            ElementIndex = -1;
        }
    }

    private void moveToPositionAligned(int position)
    {
        if(position <= RefPosition)
            return;

        int posDifference = position - RefPosition;
        int remainingElement = BaseLength - ElementIndex - 1;

        if(posDifference > remainingElement)
        {
            ReadIndex = BaseLength;
            return;
        }

        RefPosition += posDifference;
        ReadIndex += posDifference;
        ElementIndex += posDifference;
    }

    public boolean isValid()
    {
        if(Elements == null)
            return ReadIndex < BaseLength;

        return CigarIndex >= 0 && CigarIndex < Elements.size() && ReadIndex < BaseLength;
    }

    public boolean isDelete()
    {
        if(!isValid() || Elements == null)
            return false;

        return Elements.get(CigarIndex).getOperator() == D;
    }

    public boolean isInsert()
    {
        if(!isValid() || Elements == null)
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
        CigarElement element = Elements != null && isValid() ? Elements.get(CigarIndex) : null;

        return format("cigar(%d of %d) element(%d of %s) readIndex(%d) refPos(%d)",
                CigarIndex, Elements != null ? Elements.size() : 1, ElementIndex, element, ReadIndex, RefPosition);
    }
}
