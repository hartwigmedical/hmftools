package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.Math.min;
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
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class IndelConsensusReadsv1
{
    private final BaseBuilder mBaseBuilder;

    public IndelConsensusReadsv1(final BaseBuilder baseBuilder)
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

        // purge any reads with 2+ indel differences from the consensus
        final List<SAMRecord> consensusReads = purgeExcessIndelReads(reads);

        List<ReadParseState> readStates = consensusReads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int readCount = consensusReads.size();
        int minRequiredReads = (readCount % 2) == 0 ? readCount / 2 : readCount / 2 + 1;

        int baseIndex = 0;
        boolean addedTrailingSoftClip = false;

        while(!readStates.isEmpty())
        {
            if(addedTrailingSoftClip)
            {
                readStates.stream().filter(x -> x.elementType() != S).forEach(x -> x.markExhausted());
            }

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

            int activeReadCount = (int)readStates.stream().filter(x -> !x.exhausted()).count();

            // exit if there are now less than half the reads still active
            if(activeReadCount < minRequiredReads && alignedOrSoftClip(selectedElement.Operator))
                break;

            if(!deleteOrSplit(selectedElement.Operator))
            {
                if(consensusState.IsForward)
                    baseIndex += selectedElement.Length;
            }

            if(consensusState.CigarElements.size() > 1 && selectedElement.Operator == S)
                addedTrailingSoftClip = true;
        }

        resetBoundaries(consensusState);

        consensusState.setOutcome(INDEL_MISMATCH);
    }

    private void resetBoundaries(final ConsensusState consensusState)
    {
        int positionLength = consensusState.CigarElements.stream()
                .filter(x -> x.getOperator().consumesReferenceBases()).mapToInt(x -> x.getLength()).sum();

        if(consensusState.IsForward)
        {
            consensusState.MaxAlignedPosEnd = consensusState.MinAlignedPosStart + positionLength - 1;
            CigarElement lastElement = consensusState.CigarElements.get(consensusState.CigarElements.size() - 1);

            if(lastElement.getOperator() == S)
                consensusState.MaxUnclippedPosEnd = consensusState.MaxAlignedPosEnd + lastElement.getLength();
            else
                consensusState.MaxUnclippedPosEnd = consensusState.MaxAlignedPosEnd;
        }
        else
        {
            consensusState.MinAlignedPosStart = consensusState.MaxAlignedPosEnd - positionLength + 1;

            CigarElement firstElement = consensusState.CigarElements.get(0);

            if(firstElement.getOperator() == S)
                consensusState.MinUnclippedPosStart = consensusState.MinAlignedPosStart - firstElement.getLength();
            else
                consensusState.MinUnclippedPosStart = consensusState.MinAlignedPosStart;
        }
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

                // if(operatorsDiffer(read.elementType(), selectedElement.Operator))
                if(read.elementType() != selectedElement.Operator)
                {
                    // when aligned (M) is selected:
                    // - insert - skip past the insert's bases
                    // - delete - move along in step but cannot use base
                    // - ignore read with soft-clipped bases since they are likely misaligned

                    // when insert (I) is selected:
                    // - aligned or delete - pause the index and don't use the base

                    if(selectedElement.Operator == M)
                    {
                        if(deleteOrSplit(read.elementType()) || read.elementType() == S)
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

    public static boolean alignedOrSoftClip(final CigarOperator operator)
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

    private List<SAMRecord> purgeExcessIndelReads(final List<SAMRecord> reads)
    {
        Map<Integer,Integer> indelCounts = Maps.newHashMap();

        int maxIndelCount = 0;
        int minIndelCount = -1;
        int maxFreqIndelCount = 0;
        int maxFreq = 0;
        for(SAMRecord read : reads)
        {
            int indelCount = indelCount(read);
            maxIndelCount = max(maxIndelCount, indelCount);

            if(minIndelCount == -1)
                minIndelCount = indelCount;
            else
                minIndelCount = min(minIndelCount, indelCount);

            Integer count = indelCounts.get(indelCount);
            int freq = count != null ? count : 1;
            indelCounts.put(indelCount, freq);

            if(freq > maxFreq)
            {
                maxFreq = freq;
                maxFreqIndelCount = indelCount;
            }
        }

        if(indelCounts.size() == 1 || maxFreq == 1)
            return reads;

        if(maxIndelCount <= maxFreqIndelCount + 1 && minIndelCount >= maxFreqIndelCount -1)
            return reads;

        List<SAMRecord> consensusReads = Lists.newArrayListWithCapacity(reads.size());

        for(SAMRecord read : reads)
        {
            int indelCount = indelCount(read);
            if(indelCount >= maxFreqIndelCount - 1 && indelCount <= maxFreqIndelCount + 1)
                consensusReads.add(read);
        }

        return consensusReads;
    }

    private static int indelCount(final SAMRecord read)
    {
        return (int)read.getCigar().getCigarElements().stream().filter(x -> x.getOperator().isIndel()).count();
    }

    private ElementTypeCount findNextElement(final List<ReadParseState> readStates)
    {
        List<ElementTypeCount> counts = Lists.newArrayList();

        boolean allLastAlignOrSplit = true;
        boolean hasLastAligned = true;

        for(ReadParseState read : readStates)
        {
            if(read.exhausted())
                continue;

            if(read.isLastAlignOrSplit())
            {
                hasLastAligned |= read.elementType() == M;
            }
            else
            {
                allLastAlignOrSplit = false;
            }

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
        boolean allInitial = readStates.stream().allMatch(x -> x.elementIndex() == 0);
        boolean favourAligned = allInitial || allLastAlignOrSplit || hasLastAligned;

        for(ElementTypeCount elementCount : counts)
        {
            if(maxElement != null)
            {
                if(favourAligned && maxElement.Operator == S && elementCount.Operator == M)
                {
                    // favour anything over soft-clips, even Ds and Is since they'll be followed by Ms
                    maxElement = elementCount;
                }
                else if(favourAligned && elementCount.Operator == S && maxElement.Operator == M)
                {
                    // keep non-soft-clipped selection
                }
                if(elementCount.Count > maxElement.Count)
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
        int cigarCount = refCigar.getCigarElements().size();

        CigarElement firstAlignment = null;
        CigarElement lastAlignment = null;

        boolean hasFirstAlignment = refCigar.getCigarElements().get(0).getOperator() == M
                || (refCigar.getCigarElements().get(0).getOperator() == S && refCigar.getCigarElements().get(1).getOperator() == M);

        boolean hasLastAlignment = refCigar.getCigarElements().get(cigarCount - 1).getOperator() == M
                || (refCigar.getCigarElements().get(cigarCount - 1).getOperator() == S && refCigar.getCigarElements().get(cigarCount - 2).getOperator() == M);

        // indel CIGARs can occasionally begin with an 'I', so cannot assume

        for(int i = 0; i < refCigar.getCigarElements().size(); ++i)
        {
            CigarElement element = refCigar.getCigarElements().get(i);

            if(element.getOperator() == S)
                continue;
            else if(hasFirstAlignment && firstAlignment == null && element.getOperator() == M)
                firstAlignment = element;
            else if(hasLastAlignment && i >= cigarCount - 2 && element.getOperator() == M)
                lastAlignment = element;
            else
                consensusState.CigarElements.add(element);
        }

        // lengthen the outer alignments as required
        if(firstAlignment != null)
        {
            int firstAlignmentLength =
                    firstAlignment.getLength() + max(initialRead.getAlignmentStart() - consensusState.MinAlignedPosStart, 0);

            consensusState.CigarElements.add(0, new CigarElement(firstAlignmentLength, M));

            if(leftSoftClipBases > 0)
                consensusState.CigarElements.add(0, new CigarElement(leftSoftClipBases, S));
        }

        if(lastAlignment != null)
        {
            int lastAlignmentLength = lastAlignment.getLength() + max(consensusState.MaxAlignedPosEnd - initialRead.getAlignmentEnd(), 0);

            consensusState.CigarElements.add(new CigarElement(lastAlignmentLength, M));

            if(rightSoftClipBases > 0)
                consensusState.CigarElements.add(new CigarElement(rightSoftClipBases, S));
        }
    }
}
