package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarStringFromElements;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.sage.common.SimpleVariant.isLongInsert;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ReadCigarInfo
{
    public List<CigarElement> Cigar;

    public final int ReadAlignmentStart;

    // positions may be an unclipped positions if the required indices fall in soft-clips
    public final int FlankPositionStart;  // may be an unclipped positions if the required indices fall in soft-clips
    public final int FlankPositionEnd;
    public final int CorePositionStart;  // may be an unclipped positions if the required indices fall in soft-clips
    public final int CorePositionEnd;

    // calculated for the scenario where an indel in the flanks pushes out the alignment beyond the standard flank length
    public final int FlankIndexStart;
    public final int FlankIndexEnd;

    public ReadCigarInfo(
            final int readAlignmentStart, final List<CigarElement> cigar, final int flankPositionStart, final int flankPositionEnd,
            final int corePositionStart, final int corePositionEnd, final int flankIndexStart, final int flankIndexEnd)
    {
        ReadAlignmentStart = readAlignmentStart;
        Cigar = cigar;
        FlankPositionStart = flankPositionStart;
        FlankPositionEnd = flankPositionEnd;
        CorePositionStart = corePositionStart;
        CorePositionEnd = corePositionEnd;
        FlankIndexStart = flankIndexStart;
        FlankIndexEnd = flankIndexEnd;
    }

    public boolean isValid()
    {
        return FlankIndexStart >= 0 && FlankIndexEnd > FlankIndexStart
            && CorePositionStart < CorePositionEnd && CorePositionStart > FlankPositionStart && CorePositionEnd < FlankPositionEnd;
    }

    public String toString()
    {
        return format("%s flankPos(%d-%d) corePos(%d-%d) flankIndex(%d-%d)",
                cigarStringFromElements(Cigar), FlankPositionStart, FlankPositionEnd, CorePositionStart, CorePositionEnd,
                FlankIndexStart, FlankIndexEnd);
    }

    public static ReadCigarInfo buildReadCigar(
            int readStart, final List<CigarElement> cigarElements, int leftFlankIndex, int leftCoreIndex, int rightCoreIndex,
            int rightFlankIndex)
    {
        // find the read index of the variant, and then build out the read cigar from there to the flanks, noting the flank
        // and core positions along the way
        List<CigarElement> cigar = Lists.newArrayList();

        int readIndex = 0;
        int refPosition = readStart;
        int flankPosStart = 0;
        int flankPosEnd = 0;
        int corePosStart = 0;
        int corePosEnd = 0;

        int finalIndexStart = leftFlankIndex;
        int finalIndexEnd = rightFlankIndex;

        for(CigarElement element : cigarElements)
        {
            if(readIndex == 0 && element.getOperator() == S)
            {
                // set to unclipped ref position so the alignment can capture corresponding ref bases
                refPosition -= element.getLength();
            }

            int elementEndIndex = readIndex + element.getLength() - 1;
            boolean isReadBases = element.getOperator().consumesReadBases();

            if(isReadBases && readIndex <= leftFlankIndex && elementEndIndex >= rightFlankIndex)
            {
                // this element straddles all required indices, must by definition be an aligned section
                int elementStart = max(readIndex, leftFlankIndex);
                int elementEnd = min(elementEndIndex, rightFlankIndex);

                // cap this new element to what is required to reach the end of the core (ie the index end) unless it falls in an indel or split
                int elementLength = elementEnd - elementStart + 1;
                cigar.add(new CigarElement(elementLength, element.getOperator()));

                flankPosStart = refPosition + (leftFlankIndex - readIndex);
                flankPosEnd = refPosition + (rightFlankIndex - readIndex);
                corePosStart = refPosition + (leftCoreIndex - readIndex);
                corePosEnd = refPosition + (rightCoreIndex - readIndex);
                break;
            }

            boolean addedElement = false;

            if(elementEndIndex >= leftFlankIndex)
            {
                int readEndIndex = isReadBases ? elementEndIndex : readIndex;

                // first check setting the core positions
                if(corePosStart == 0 && positionWithin(leftCoreIndex, readIndex, readEndIndex))
                {
                    // adjust position back if in an indel
                    if(element.getOperator() == I)
                        corePosStart = refPosition - 1;
                    else if(element.getOperator() == D || element.getOperator() == N)
                        corePosStart = refPosition + element.getLength(); // read index points to the next non-skipped base
                    else
                        corePosStart = refPosition + (leftCoreIndex - readIndex);
                }

                if(corePosEnd == 0 && positionWithin(rightCoreIndex, readIndex, readEndIndex))
                {
                    // adjust position forward if in an indel
                    if(element.getOperator() == I)
                    {
                        corePosEnd = refPosition; // already pointing at the next M (aligned) base
                    }
                    else if(element.getOperator() == D || element.getOperator() == N)
                    {
                        corePosEnd = refPosition + element.getLength();
                    }
                    else
                    {
                        corePosEnd = refPosition + (rightCoreIndex - readIndex);
                    }
                }

                if(isReadBases)
                {
                    if(flankPosStart == 0)
                    {
                        // covers the lower required index only - check if this falls in an indel
                        int elementStart = max(readIndex, leftFlankIndex);
                        int elementEnd = elementEndIndex;

                        // cap this new element to what is required to reach the end of the core (ie the index end) unless it falls in an indel or split
                        int elementLength = element.getOperator().isIndel() || element.getOperator() == N ?
                                element.getLength() : elementEnd - elementStart + 1;

                        cigar.add(new CigarElement(elementLength, element.getOperator()));

                        if(element.getOperator().isIndel())
                        {
                            // handles an insert that pushes the alignment out - always take the prior alignment base and reduce index start
                            // eg looking to find the alignment boundary for index start 12 for 10M5I1350M, alignment start == 100
                            // so at the insert element, read index = 10, ref pos = 110 (pointing at next ref base)
                            cigar.remove(0);
                            int extraIndexStart = leftFlankIndex - readIndex + 1;
                            cigar.add(0, new CigarElement(1, M));

                            finalIndexStart -= extraIndexStart;
                            flankPosStart = max(refPosition - 1, readStart);
                        }
                        else
                        {
                            flankPosStart = refPosition + (leftFlankIndex - readIndex);
                        }

                        addedElement = true;
                    }
                    else if(flankPosEnd == 0 && elementEndIndex >= rightFlankIndex)
                    {
                        int elementStart = readIndex;

                        if(element.getOperator().isIndel() || element.getOperator() == N)
                        {
                            cigar.add(element);
                            cigar.add(new CigarElement(1, M));

                            if(element.getOperator() == I)
                            {
                                flankPosEnd = refPosition; // already pointing at the next M (aligned) base

                                // similar extension to the above
                                int extraIndexEnd = elementEndIndex + 1 - rightFlankIndex;
                                finalIndexEnd += extraIndexEnd;
                            }
                            else
                            {
                                flankPosEnd = refPosition + element.getLength();
                            }
                        }
                        else
                        {
                            int elementEnd = min(elementEndIndex, rightFlankIndex);
                            int elementLength = elementEnd - elementStart + 1;
                            cigar.add(new CigarElement(elementLength, element.getOperator()));

                            flankPosEnd = refPosition + (rightFlankIndex - readIndex);
                        }

                        break;
                    }
                }
            }

            if(!addedElement)
            {
                // add this element falling in between the required indices
                cigar.add(element);
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
                refPosition += element.getLength();

            if(readIndex > rightFlankIndex)
                break;
        }

        return new ReadCigarInfo(readStart, cigar, flankPosStart, flankPosEnd, corePosStart, corePosEnd, finalIndexStart, finalIndexEnd);
    }
}
