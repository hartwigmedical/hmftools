package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarStringFromElements;

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
    public final int UnclippedStart;  // unclipped positions, ie possibly extending past the read's alignment boundaries
    public final int UnclippedEnd;

    // calculated for the scenario where an indel in the flanks pushes out the alignment beyond the standard flank length
    public final int FlankIndexStart;
    public final int FlankIndexEnd;

    public ReadCigarInfo(
            final List<CigarElement> cigar, final int unclippedStart, final int unclippedEnd,
            final int flankIndexStart, final int flankIndexEnd)
    {
        Cigar = cigar;
        UnclippedStart = unclippedStart;
        UnclippedEnd = unclippedEnd;
        FlankIndexStart = flankIndexStart;
        FlankIndexEnd = flankIndexEnd;
    }

    public String toString()
    {
        return format("%s align(%d-%d) flankIndex(%d-%d)",
                cigarStringFromElements(Cigar), UnclippedStart, UnclippedEnd, FlankIndexStart, FlankIndexEnd);
    }

    public static ReadCigarInfo buildReadCigar(final SAMRecord read, int indexStart, int indexEnd)
    {
        // builds a cigar around the specified index boundaries and calculates the corresponding alignment positions
        List<CigarElement> cigar = Lists.newArrayList();

        int readIndex = 0;
        int refPosition = read.getAlignmentStart();
        int unclippedPosStart = 0;
        int unclippedPosEnd = 0;

        int finalIndexStart = indexStart;
        int finalIndexEnd = indexEnd;

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            if(readIndex == 0 && element.getOperator() == S)
            {
                // set to unclipped ref position so the alignment can capture corresponding ref bases
                refPosition -= element.getLength();
            }

            int elementEndIndex = readIndex + element.getLength() - 1;

            if(elementEndIndex >= indexStart)
            {
                int elementStart = max(readIndex, indexStart);
                int elementEnd = min(elementEndIndex, indexEnd);

                // cap this new element to what is required to reach the end of the core (ie the index end) unless it falls in an indel or split
                int elementLength = element.getOperator().isIndel() || element.getOperator() == N ?
                        element.getLength() : elementEnd - elementStart + 1;

                cigar.add(new CigarElement(elementLength, element.getOperator()));

                if(unclippedPosStart == 0)
                {
                    if(element.getOperator().isIndel())
                    {
                        // handles an insert that pushes the alignment out - always take the prior alignment base and reduce index start
                        // eg looking to find the alignment boundary for index start 12 for 10M5I1350M, alignment start == 100
                        // so at the insert element, read index = 10, ref pos = 110 (pointing at next ref base)
                        int extraIndexStart = indexStart - readIndex + 1;
                        cigar.add(0, new CigarElement(1, M));

                        finalIndexStart -= extraIndexStart;
                        unclippedPosStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    /*
                    else if(element.getOperator() == D)
                    {
                        unclippedPosStart = max(refPosition - 1, read.getAlignmentStart());
                    }
                    */
                    else
                    {
                        unclippedPosStart = refPosition + (indexStart - readIndex);
                    }
                }

                if(unclippedPosEnd == 0 && elementEndIndex >= indexEnd)
                {
                    if(element.getOperator() == I)
                    {
                        // similar extension to the above
                        int extraIndexEnd = elementEndIndex + 1 - indexEnd;
                        cigar.add(new CigarElement(1, M));

                        finalIndexEnd += extraIndexEnd;

                        unclippedPosEnd = refPosition; // already pointing at the next M (aligned) base
                    }
                    else if(element.getOperator() == D)
                    {
                        unclippedPosEnd = refPosition + element.getLength();
                    }
                    else
                    {
                        unclippedPosEnd = refPosition + (indexEnd - readIndex); // unclipped so not bound by the read's alignment
                    }
                }
            }

            if(element.getOperator().consumesReadBases())
                readIndex += element.getLength();

            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
                refPosition += element.getLength();

            if(readIndex > indexEnd)
                break;
        }

        return new ReadCigarInfo(cigar, unclippedPosStart, unclippedPosEnd, finalIndexStart, finalIndexEnd);
    }
}
