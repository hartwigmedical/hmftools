package com.hartwig.hmftools.common.bam;

import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class ReadCigarState
{
    public int RefPosition;
    public int ReadIndex;
    public int CigarIndex;
    public int ElementIndex;
    public CigarElement Element;

    private boolean mValid;

    public ReadCigarState(
            final int refPosition, final int readIndex, final CigarElement element, final int cigarIndex, final int elementIndex)
    {
        RefPosition = refPosition;
        ReadIndex = readIndex;
        Element = element;
        CigarIndex = cigarIndex;
        ElementIndex = elementIndex;
        mValid = true;
    }

    public ReadCigarState(final ReadCigarState other)
    {
        this(other.RefPosition, other.ReadIndex, other.Element, other.CigarIndex, other.ElementIndex);
    }

    public boolean isValid() { return mValid; }
    public void resetValid() { mValid = true; }
    public void setInvalid() { mValid = false; }

    public static ReadCigarState initialise(final int refPosition, final List<CigarElement> cigarElements)
    {
        if(cigarElements.get(0).getOperator() == S)
        {
            // move past initial soft-clip
            return new ReadCigarState(refPosition, cigarElements.get(0).getLength(), cigarElements.get(1), 1, 0);
        }
        else
        {
            return new ReadCigarState(refPosition, 0, cigarElements.get(0), 0, 0);
        }
    }

    public CigarOperator operator()
    {
        return Element.getOperator();
    }

    public String toString()
    {
        return format("ref(%d) read(%d) cigar(%d: %d/%d%s)%s",
                RefPosition, ReadIndex, CigarIndex, ElementIndex + 1, Element.getLength(), Element.getOperator(),
                mValid ? "" : " invalid");
    }

    // navigation methods
    public static void moveState(final ReadCigarState state, final List<CigarElement> cigarElements, boolean moveUp)
    {
        // exit on any soft-clipped bases or once at the end of the read's cigar
        CigarOperator lastOperator = state.operator();

        if(moveUp)
        {
            ++state.ElementIndex;

            if(state.ElementIndex >= state.Element.getLength())
            {
                ++state.CigarIndex;
                if(state.CigarIndex >= cigarElements.size())
                {
                    state.setInvalid();
                    return;
                }

                state.Element = cigarElements.get(state.CigarIndex);
                state.ElementIndex = 0;
            }

            if(lastOperator != S && state.operator().consumesReferenceBases()) // since soft-clip will have remained at the read start
                ++state.RefPosition;

            if(state.operator().consumesReadBases())
                ++state.ReadIndex;
        }
        else
        {
            --state.ElementIndex;

            if(state.ElementIndex < 0)
            {
                --state.CigarIndex;
                if(state.CigarIndex < 0)
                {
                    state.setInvalid();
                    return;
                }

                state.Element = cigarElements.get(state.CigarIndex);
                state.ElementIndex = state.Element.getLength() - 1;
            }

            if(lastOperator != S && state.operator().consumesReferenceBases())
                --state.RefPosition;

            if(state.operator().consumesReadBases())
                --state.ReadIndex;
        }
    }

    public static void moveToRefPosition(final ReadCigarState state, final List<CigarElement> cigarElements, final int targetRefPosition)
    {
        boolean moveUp = state.RefPosition < targetRefPosition;

        while(state.RefPosition != targetRefPosition)
        {
            moveState(state, cigarElements, moveUp);

            if(!state.isValid())
                break;
        }

        if(state.RefPosition != targetRefPosition)
            state.setInvalid();
    }

    public static void moveToIndex(final ReadCigarState state, final List<CigarElement> cigarElements, final int targetIndex)
    {
        while(state.ReadIndex != targetIndex)
        {
            moveState(state, cigarElements, state.ReadIndex < targetIndex);

            if(!state.isValid())
                break;
        }

        if(state.ReadIndex != targetIndex)
            state.setInvalid();
    }
}
