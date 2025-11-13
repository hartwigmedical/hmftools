package com.hartwig.hmftools.sage.seqtech;

import static com.hartwig.hmftools.common.bam.ReadCigarState.moveToIndex;
import static com.hartwig.hmftools.common.bam.ReadCigarState.moveToRefPosition;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ReadCigarState;
import com.hartwig.hmftools.sage.common.ReadCigarInfo;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class UltimaCoreExtender
{
    public static ReadCigarInfo extendUltimaCore(
            final byte[] readBases, final RefSequence refSequence, final int readAlignmentStart,
            final List<CigarElement> cigarElements, final ReadCigarInfo readCigarInfo, final int flankSize, final boolean inAppendMode)
    {
        int corePositionStart = readCigarInfo.CorePositionStart;
        int corePositionEnd = readCigarInfo.CorePositionEnd;
        int flankPositionStart = readCigarInfo.FlankPositionStart;
        int flankIndexStart = readCigarInfo.FlankIndexStart;
        int flankPositionEnd = readCigarInfo.FlankPositionEnd;
        int flankIndexEnd = readCigarInfo.FlankIndexEnd;

        // move to the start of the core
        ReadCigarState lowerState = ReadCigarState.initialise(readAlignmentStart, cigarElements);
        ReadCigarState upperState = new ReadCigarState(lowerState);

        int readAlignmentEnd = readAlignmentStart + cigarElements.stream()
                .filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum() - 1;

        boolean extendedStart = false;
        boolean extendedEnd = false;

        boolean lowerInSoftClip = readCigarInfo.FlankPositionStart < readAlignmentStart;
        boolean upperInSoftClip = readCigarInfo.FlankPositionEnd > readAlignmentEnd;

        if(!lowerInSoftClip)
        {
            moveToRefPosition(lowerState, cigarElements, readCigarInfo.CorePositionStart);

            if(lowerState.isValid())
                extendedStart = findCoreExtension(readBases, refSequence, cigarElements, lowerState, false);

            if(extendedStart)
            {
                corePositionStart = lowerState.RefPosition;

                findFlankIndex(lowerState, cigarElements, flankSize, false);

                if(!lowerState.isValid() && inAppendMode)
                {
                    // move to the first aligned base, ie up from the soft-clip or end of the read
                    lowerState.resetValid();
                    moveState(lowerState, cigarElements, true, true);
                }

                flankPositionStart = lowerState.RefPosition;
                flankIndexStart = lowerState.ReadIndex;
            }
            else
            {
                // keep the existing flank position
                lowerState = ReadCigarState.initialise(readAlignmentStart, cigarElements);
                moveToRefPosition(lowerState, cigarElements, readCigarInfo.FlankPositionStart);
            }
        }
        else
        {
            // move to the existing flank start
            moveToIndex(lowerState, cigarElements, readCigarInfo.FlankIndexStart);
        }

        if(!upperInSoftClip)
        {
            moveToRefPosition(upperState, cigarElements, readCigarInfo.CorePositionEnd);

            if(upperState.isValid())
                extendedEnd = findCoreExtension(readBases, refSequence, cigarElements, upperState, true);

            if(extendedEnd)
            {
                corePositionEnd = upperState.RefPosition;

                findFlankIndex(upperState, cigarElements, flankSize, true);

                if(!upperState.isValid() && inAppendMode)
                {
                    upperState.resetValid();
                    moveState(upperState, cigarElements, false, true);
                }

                flankPositionEnd = upperState.RefPosition;
                flankIndexEnd = upperState.ReadIndex;
            }
            else
            {
                moveToRefPosition(upperState, cigarElements, readCigarInfo.FlankPositionEnd);
            }
        }
        else
        {
            moveToIndex(upperState, cigarElements, readCigarInfo.FlankIndexEnd);
        }

        if(!extendedStart && !extendedEnd)
            return readCigarInfo;

        if(!lowerState.isValid() || !upperState.isValid())
            return null;

        // build a new cigar between the 2 flank boundaries
        List<CigarElement> newCigarElements = buildReadCigar(new ReadCigarState(lowerState), cigarElements, upperState.ReadIndex);

        return new ReadCigarInfo(
                readAlignmentStart, newCigarElements, flankPositionStart, flankPositionEnd, corePositionStart, corePositionEnd,
                flankIndexStart, flankIndexEnd);
    }

    private static void moveState(final ReadCigarState state, final List<CigarElement> cigarElements, boolean moveUp, boolean skipDeletes)
    {
        ReadCigarState.moveState(state, cigarElements, moveUp);

        while(skipDeletes && state.operator() == D)
        {
            ReadCigarState.moveState(state, cigarElements, moveUp);
        }

        if(state.operator() == S) // for this routine, the adjusted core and flanks cannot be within a soft-clip
            state.setInvalid();
    }

    private static boolean findCoreExtension(
            final byte[] readBases, final RefSequence refSequence, final List<CigarElement> cigarElements,
            final ReadCigarState state, boolean searchUp)
    {
        // extend the core while any homopolymer exists in either ref or base, then exit when the ref matches the base at an M element
        byte readBase = readBases[state.ReadIndex];
        int refPosition = state.RefPosition;
        byte refBase = refSequence.base(refPosition);

        boolean requiredExtension = false;

        while(true)
        {
            moveState(state, cigarElements, searchUp, true);
            refPosition += searchUp ? 1 : -1;

            if(!state.isValid())
                break;

            byte nextReadBase = readBases[state.ReadIndex];

            if(!refSequence.containsPosition(refPosition))
            {
                state.setInvalid();
                return false;
            }

            byte nextRefBase = refSequence.base(refPosition);

            boolean hasHomopolymer = (readBase == nextReadBase) || (refBase == nextRefBase);
            boolean nextBasesMatch = nextReadBase == nextRefBase;

            if(!requiredExtension && !hasHomopolymer && nextBasesMatch) // immediate exit if no extension was required
                return false;

            if(hasHomopolymer)
            {
                requiredExtension = true;
            }
            else
            {
                if(nextBasesMatch && (state.operator() == CigarOperator.M || state.operator() == CigarOperator.S))
                    break;
            }

            readBase = nextReadBase;
            refBase = nextRefBase;
        }

        return requiredExtension;
    }

    private static void findFlankIndex(
            final ReadCigarState state, final List<CigarElement> cigarElements, final int requiredFlankLength, boolean moveUp)
    {
        if(state.operator() == S)
            return;

        int flankLength = 0;

        while(flankLength < requiredFlankLength || state.Element.getOperator() == I)
        {
            moveState(state, cigarElements, moveUp, true);
            ++flankLength;

            if(!state.isValid())
                break;
        }
    }

    private static List<CigarElement> buildReadCigar(
            final ReadCigarState state, final List<CigarElement> cigarElements, final int flankIndexEnd)
    {
        List<CigarElement> elements = Lists.newArrayList();
        int cigarLength = 1;
        CigarOperator cigarOperator = state.operator();

        while(state.ReadIndex < flankIndexEnd)
        {
            moveState(state, cigarElements, true, false);

            if(state.operator() == cigarOperator)
            {
                ++cigarLength;
            }
            else
            {
                elements.add(new CigarElement(cigarLength, cigarOperator));
                cigarLength = 1;
                cigarOperator = state.operator();
            }
        }

        elements.add(new CigarElement(cigarLength, cigarOperator));
        return elements;
    }
}
