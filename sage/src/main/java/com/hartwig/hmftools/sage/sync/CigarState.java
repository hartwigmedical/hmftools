package com.hartwig.hmftools.sage.sync;

import static java.lang.Math.min;
import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class CigarState
{
    public final List<CigarElement> Elements;
    public final int ElementCount;
    public int ElementIndex; // index of the current cigar element in the list
    public CigarElement CurrentElement;
    public int ElementBaseIndex; // index of the current base within the current element

    public CigarState(List<CigarElement> elements)
    {
        Elements = elements;
        ElementCount = elements.size();
        ElementIndex = 0;
        CurrentElement = Elements.get(ElementIndex);
        ElementBaseIndex = 0;
    }

    public boolean exhausted() { return CurrentElement == null; }
    public int remainingElementBases() { return CurrentElement.getLength() - ElementBaseIndex; }
    public CigarOperator currentOperator() { return CurrentElement != null ? CurrentElement.getOperator() : null; }

    public int softClipStart() { return Elements.get(0).getOperator() == S ? Elements.get(0).getLength() : 0; }
    public int softClipEnd() { return Elements.get(ElementCount - 1).getOperator() == S ? Elements.get(ElementCount - 1).getLength() : 0; }

    public boolean moveNext()
    {
        // returns true if changes element or ends
        if(exhausted())
            return false;

        if(ElementBaseIndex < CurrentElement.getLength() - 1)
        {
            ++ElementBaseIndex;
            return false;
        }

        if(ElementIndex == ElementCount - 1)
        {
            CurrentElement = null;
        }
        else
        {
            ++ElementIndex;
            CurrentElement = Elements.get(ElementIndex);
            ElementBaseIndex = 0;
        }

        return true;
    }

    public String toString()
    {
        if(exhausted())
        {
            return format("cigar(%d/%d) exhausted", ElementIndex, ElementCount);
        }
        else
        {
            return format("cigar(%d/%d) element(%s %d/%d)",
                    ElementIndex, ElementCount, CurrentElement.getOperator(), ElementBaseIndex, CurrentElement.getLength());
        }
    }

    public void moveToPosition(int effectiveStart, int requiredStart)
    {
        int positionDiff = requiredStart - effectiveStart;

        while(positionDiff > 0)
        {
            int minBases = min(CurrentElement.getLength(), positionDiff);

            if(CurrentElement.getOperator().consumesReferenceBases() || CurrentElement.getOperator() == S)
                positionDiff -= minBases;

            if(minBases < CurrentElement.getLength())
            {
                ElementBaseIndex += minBases;
                break;
            }
            else
            {
                ++ElementIndex;
                CurrentElement = Elements.get(ElementIndex);
                ElementBaseIndex = 0;
            }

            if(positionDiff <= 0)
                break;
        }
    }
}
