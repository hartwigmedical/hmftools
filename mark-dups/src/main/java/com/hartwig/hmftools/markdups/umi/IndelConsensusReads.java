package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;

import static com.hartwig.hmftools.markdups.umi.BaseBuilder.NO_BASE;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MATCH;

import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;
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
            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(INDEL_MATCH);
            buildIndelCigar(consensusState,reads.get(0));
            return;
        }

        List<ReadParseState> readStates = reads.stream().map(x -> new ReadParseState(x, consensusState.IsForward)).collect(Collectors.toList());

        int readCount = reads.size();
        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        int baseIndex = consensusState.IsForward ? 0 : consensusState.baseLength() - 1;
        CigarOperator currentOperator = M;
        int currentCigarLength = 0;

        int exhaustedCount = 0;
        ElementTypeCount selectedElement = null;

        while(exhaustedCount < readCount)
        {
            // work out current element type
            if(selectedElement == null)
                selectedElement = findNextElement(readStates);

            // simplest scenario is where all reads agree about this next element
            if(selectedElement.Count == readCount)
            {
                addMatchingElement(consensusState, readStates, selectedElement, baseIndex);

                for(ReadParseState read : readStates)
                {
                    if(read.exhausted())
                        ++exhaustedCount;
                }

                selectedElement = null;

                if(consensusState.IsForward)
                    baseIndex += selectedElement.Length;
                else
                    baseIndex -= selectedElement.Length;

                continue;
            }

            //


            for(ReadParseState read : readStates)
            {

            }




        }

    }

    private void addMatchingElement(
            final ConsensusState consensusState, final List<ReadParseState> readStates, final ElementTypeCount selectedElement, int baseIndex)
    {
        int readCount = readStates.size();
        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        consensusState.addCigarElement(selectedElement.Length, selectedElement.Operator);

        for(int i = 0; i < selectedElement.Length; ++i)
        {
            boolean hasMismatch = false;
            int maxQual = 0;

            locationBases[0] = readStates.get(0).currentBase();
            locationQuals[0] = readStates.get(0).currentBaseQual();
            byte firstBase = locationBases[0];

            for(int r = 1; r < readCount; ++r)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                ReadParseState read = readStates.get(r);

                locationBases[r] = read.currentBase();
                locationQuals[r] = read.currentBaseQual();

                hasMismatch |= locationBases[r] != firstBase;

                maxQual = max(locationQuals[r], maxQual);
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

        }
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
                    if(readElement.getOperator() != M && readElement.getOperator() != S)
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
            }
        }

        return true;
    }

    private ElementTypeCount findNextElement(final List<ReadParseState> readStates)
    {
        List<ElementTypeCount> counts = Lists.newArrayList();

        for(ReadParseState read : readStates)
        {
            int index = 0;
            while(index < counts.size())
            {
                ElementTypeCount elementCount = counts.get(index);

                if(elementCount.Operator == read.currentElementType() && elementCount.Length == read.currentElementLength())
                {
                    ++elementCount.Count;
                    break;
                }

                ++index;
            }

            if(index >= counts.size())
            {
                counts.add(new ElementTypeCount(read.currentElementType(), read.currentElementLength()));
            }
        }

        ElementTypeCount maxElement = null;

        for(ElementTypeCount elementCount : counts)
        {
            if(maxElement == null)
            {
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
            Count = 0;
        }
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
