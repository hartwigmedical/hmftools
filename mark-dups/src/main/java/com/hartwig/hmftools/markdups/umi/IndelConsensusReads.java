package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.umi.BaseBuilder.NO_BASE;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_FAIL;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MATCH;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MISMATCH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class IndelConsensusReads
{
    private final BaseBuilder mBaseBuilder;

    public IndelConsensusReads(final BaseBuilder baseBuilder)
    {
        mBaseBuilder = baseBuilder;
    }

    public void buildIndelComponents(final List<SAMRecord> reads, final ConsensusState consensusState)
    {
        if(haveConsistentCigars(reads))
        {
            int baseLength = reads.get(0).getReadBases().length;
            consensusState.setBaseLength(baseLength);

            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(INDEL_MATCH);
            buildIndelCigar(consensusState,reads.get(0));
            return;
        }

        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int readCount = reads.size();

        int baseIndex = 0;

        int exhaustedCount = 0;

        while(exhaustedCount < readCount)
        {
            // work out current element type
            ElementTypeCount selectedElement = findNextElement(readStates);

            if(!deleteOrSplit(selectedElement.Operator))
            {
                consensusState.expandBaseLength(selectedElement.Length);

                if(!consensusState.IsForward)
                    baseIndex = selectedElement.Length - 1;
            }

            // simplest scenario is where all reads agree about this next element
            addElementBases(consensusState, readStates, selectedElement, baseIndex);

            if(consensusState.outcome() == INDEL_FAIL)
                return;

            for(ReadParseState read : readStates)
            {
                if(read.exhausted())
                    ++exhaustedCount;
            }

            if(!deleteOrSplit(selectedElement.Operator))
            {
                if(consensusState.IsForward)
                    baseIndex += selectedElement.Length;
            }
        }

        consensusState.setOutcome(INDEL_MISMATCH);
    }

    private void addElementBases(
            final ConsensusState consensusState, final List<ReadParseState> readStates, final ElementTypeCount selectedElement, int baseIndex)
    {
        int readCount = readStates.size();

        consensusState.addCigarElement(selectedElement.Length, selectedElement.Operator);

        if(deleteOrSplit(selectedElement.Operator))
        {
            // move past the delete element and any differing aligned bases
            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted())
                    continue;

                for(int i = 0; i < selectedElement.Length; ++i)
                {
                    if(read.elementType() == I)
                        read.skipInsert();

                    if(deleteOrSplit(read.elementType()) || alignedOrSoftClip(read.elementType()))
                    {
                        read.moveNext();
                    }
                }
            }

            return;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        for(int i = 0; i < selectedElement.Length; ++i)
        {
            boolean hasMismatch = false;
            int maxQual = 0;
            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                locationBases[r] = NO_BASE;
            }

            for(int r = 0; r < readCount; ++r)
            {
                ReadParseState read = readStates.get(r);

                if(read.exhausted())
                    continue;

                // check for element type differences:

                // first skip past any insert if the selected element is aligned
                if(selectedElement.Operator == M && read.elementType() == I)
                    read.skipInsert();

                boolean useBase = true;
                boolean moveNext = true;

                if(operatorsDiffer(read.elementType(), selectedElement.Operator))
                {
                    // when aligned (M) is selected:
                    // - insert - skip past the insert's bases
                    // - delete - move along in step but cannot use base

                    // when insert (I) is selected:
                    // - aligned or delete - pause the index and don't use the base

                    if(selectedElement.Operator == M)
                    {
                        if(deleteOrSplit(read.elementType()))
                        {
                            useBase = false;
                        }
                        else if(read.elementType() == I)
                        {
                            // handled above, implies a bug or consecutive insert
                            logMismatchFail(consensusState, r, read, selectedElement);
                            consensusState.setOutcome(INDEL_FAIL);
                            return;
                        }
                    }
                    else if(selectedElement.Operator == I)
                    {
                        if(read.elementType() == M || deleteOrSplit(read.elementType()))
                        {
                            moveNext = false;
                            useBase = false;
                        }
                    }
                }

                if(useBase)
                {
                    locationBases[r] = read.currentBase();
                    locationQuals[r] = read.currentBaseQual();

                    if(firstBase == NO_BASE)
                        firstBase = locationBases[r];
                    else
                        hasMismatch |= locationBases[r] != firstBase;

                    maxQual = max(locationQuals[r], maxQual);
                }

                if(moveNext)
                    read.moveNext();
            }

            if(!hasMismatch)
            {
                consensusState.Bases[baseIndex] = firstBase;
                consensusState.BaseQualities[baseIndex] = (byte) maxQual;
            }
            else
            {
                int basePosition = consensusState.MinUnclippedPosStart + baseIndex;

                byte[] consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(
                        locationBases, locationQuals, consensusState.Chromosome, basePosition);

                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }

            if(consensusState.IsForward)
                ++baseIndex;
            else
                --baseIndex;
        }
    }

    private void logMismatchFail(
            final ConsensusState consensusState, int readIndex, final ReadParseState read, final ElementTypeCount selectedElement)
    {
        MD_LOGGER.debug("indel mismatch fail: consensus({}:{}-{}) read({}: {}) state({}) select({})",
                consensusState.Chromosome, consensusState.MinAlignedPosStart, consensusState.MaxAlignedPosEnd,
                readIndex, read.Read.getReadName(), read.toString(), selectedElement);
    }

    private static boolean operatorsDiffer(final CigarOperator first, final CigarOperator second)
    {
        if(first == second)
            return false;

        if(alignedOrSoftClip(first) && alignedOrSoftClip(second))
            return false;

        return true;
    }

    private static boolean deleteOrSplit(final CigarOperator operator)
    {
        return operator == D || operator == N;
    }

    private static boolean alignedOrSoftClip(final CigarOperator operator)
    {
        return operator == M || operator == S;
    }

    public static boolean haveConsistentCigars(final List<SAMRecord> reads)
    {
        boolean[] hasSoftClip = new boolean[reads.size()];

        for(int i = 0; i < reads.size(); ++i)
        {
            SAMRecord read = reads.get(i);

            if(read.getCigar().isLeftClipped())
                hasSoftClip[i] = true;
        }

        SAMRecord firstRead = reads.get(0);
        int firstElementCount = firstRead.getCigar().getCigarElements().size();

        for(int c = 0; c < firstElementCount; ++c)
        {
            CigarElement firstElement = firstRead.getCigar().getCigarElements().get(c);

            if(firstElement.getOperator() == S && c == 0)
                continue;

            boolean checkLength = firstElement.getOperator() != M || (c > 1 && c < firstElementCount - 2);
            boolean isLastSoftClip = firstElement.getOperator() == S && (c == firstElementCount - 1);

            for(int i = 1; i < reads.size(); ++i)
            {
                SAMRecord read = reads.get(i);
                int cigarIndex = c;

                if(!hasSoftClip[0] && hasSoftClip[i])
                    ++cigarIndex;
                else if(hasSoftClip[0] && !hasSoftClip[i])
                    --cigarIndex;

                if(isLastSoftClip && cigarIndex >= read.getCigar().getCigarElements().size() - 1)
                    continue;

                if(!isLastSoftClip && cigarIndex >= read.getCigar().getCigarElements().size())
                    return false;

                CigarElement readElement = read.getCigar().getCigarElements().get(cigarIndex);

                if(isLastSoftClip)
                {
                    if(!alignedOrSoftClip(readElement.getOperator()))
                        return false;

                    if(cigarIndex < read.getCigar().getCigarElements().size() - 2) // must last or second last
                        return false;
                }
                else
                {
                    if(readElement.getOperator() != firstElement.getOperator())
                        return false;

                    if(checkLength && readElement.getLength() != firstElement.getLength())
                        return false;
                }

                // check for additional elements
                if(c == firstElementCount - 1 && cigarIndex < read.getCigar().getCigarElements().size() - 1)
                {
                    CigarElement nextReadElement = read.getCigar().getCigarElements().get(cigarIndex);

                    if(nextReadElement.getOperator() != S)
                        return false;
                }
            }
        }

        return true;
    }

    private ElementTypeCount findNextElement(final List<ReadParseState> readStates)
    {
        List<ElementTypeCount> counts = Lists.newArrayList();

        for(ReadParseState read : readStates)
        {
            int elementLength = read.remainingElementLength();
            int index = 0;
            while(index < counts.size())
            {
                ElementTypeCount elementCount = counts.get(index);

                if(elementCount.Operator == read.elementType() && elementCount.Length == elementLength)
                {
                    ++elementCount.Count;
                    break;
                }

                ++index;
            }

            if(index >= counts.size())
            {
                counts.add(new ElementTypeCount(read.elementType(), elementLength));
            }
        }

        ElementTypeCount maxElement = null;

        for(ElementTypeCount elementCount : counts)
        {
            if(maxElement != null)
            {
                if(elementCount.Operator == M && maxElement.Operator == S)
                {
                    // favour alignments over soft-clips
                    maxElement = elementCount;
                }
                else if(elementCount.Count > maxElement.Count)
                {
                    maxElement = elementCount;
                }
                else if(elementCount.Count == maxElement.Count
                && cigarOperatorRank(elementCount.Operator) < cigarOperatorRank(maxElement.Operator))
                {
                    maxElement = elementCount;
                }
            }
            else
            {
                maxElement = elementCount;
            }
        }

        return maxElement;
    }

    private static int cigarOperatorRank(final CigarOperator operator) { return operator.ordinal(); }

    private class ElementTypeCount
    {
        public final CigarOperator Operator;
        public final int Length;
        public int Count;

        public ElementTypeCount(final CigarOperator operator, final int length)
        {
            Operator = operator;
            Length = length;
            Count = 1;
        }

        public String toString() { return format("%d%s=%d", Length, Operator, Count); }
    }

    public void buildIndelCigar(final ConsensusState consensusState, final SAMRecord initialRead)
    {
        // build CIGAR from matched and any soft-clipped elements
        int leftSoftClipBases = consensusState.MinAlignedPosStart - consensusState.MinUnclippedPosStart;
        int rightSoftClipBases = consensusState.MaxUnclippedPosEnd - consensusState.MaxAlignedPosEnd;

        final Cigar refCigar = initialRead.getCigar();

        CigarElement firstAlignment = null;
        CigarElement lastAlignment = null;

        for(int i = 0; i < refCigar.getCigarElements().size(); ++i)
        {
            CigarElement element = refCigar.getCigarElements().get(i);

            if(element.getOperator() == S)
                continue;
            else if(firstAlignment == null && element.getOperator() == M)
                firstAlignment = element;
            else if(i >= refCigar.getCigarElements().size() - 2 && element.getOperator() == M)
                lastAlignment = element;
            else
                consensusState.CigarElements.add(element);
        }

        // lengthen the outer alignments as required
        int firstAlignmentLength = firstAlignment.getLength() + max(initialRead.getAlignmentStart() - consensusState.MinAlignedPosStart, 0);
        int lastAlignmentLength = lastAlignment.getLength() + max(consensusState.MaxAlignedPosEnd - initialRead.getAlignmentEnd(), 0);

        consensusState.CigarElements.add(0, new CigarElement(firstAlignmentLength, M));

        if(leftSoftClipBases > 0)
            consensusState.CigarElements.add(0, new CigarElement(leftSoftClipBases, S));

        consensusState.CigarElements.add(new CigarElement(lastAlignmentLength, M));

        if(rightSoftClipBases > 0)
            consensusState.CigarElements.add(new CigarElement(rightSoftClipBases, S));
    }
}
