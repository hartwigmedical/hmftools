package com.hartwig.hmftools.markdups.umi;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.ALIGNMENT_ONLY;
import static com.hartwig.hmftools.markdups.umi.ConsensusOutcome.INDEL_MATCH;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

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
        List<ReadCigarData> readCigars = Lists.newArrayListWithExpectedSize(reads.size());

        for(SAMRecord read : reads)
        {
            int readIndex = 0;
            int position = read.getAlignmentStart();

            List<CigarElementData> elementData = Lists.newArrayList();
            int cigarElementsCount = read.getCigar().getCigarElements().size();

            for(int i = 0; i < cigarElementsCount; ++i)
            {
                CigarElement element = read.getCigar().getCigarElements().get(i);

                if(element.getOperator() == S && i == 0)
                {
                    CigarElement nextElement = read.getCigar().getCigarElements().get(i + 1);
                    int unclippedPos = position - element.getLength();
                    int combinedLength = element.getLength() + nextElement.getLength();
                    elementData.add(new CigarElementData(readIndex, unclippedPos, M, combinedLength));
                    readIndex += combinedLength;
                    ++i;
                    continue;
                }
                else if(element.getOperator() == M && i == cigarElementsCount - 2)
                {
                    CigarElement nextElement = read.getCigar().getCigarElements().get(i + 1);
                    int combinedLength = element.getLength() + nextElement.getLength();
                    elementData.add(new CigarElementData(readIndex, position, M, combinedLength));
                    break;
                }

                elementData.add(new CigarElementData(readIndex, position, element.getOperator(), element.getLength()));

                if(element.getOperator() == D || element.getOperator() == N)
                {
                    position += element.getLength();
                }
                else if(element.getOperator() == M)
                {
                    readIndex += element.getLength();
                    position += element.getLength();
                }
                else if(element.getOperator() == I)
                {
                    readIndex += element.getLength() + 1;
                    position += 1;
                }
            }

            ReadCigarData readCigarData = new ReadCigarData(elementData);
            readCigars.add(readCigarData);
        }

        // check for an exact match on elements
        ReadCigarData first = readCigars.get(0);
        boolean hasMismatch = false;

        for(int i = 1; i < readCigars.size(); ++i)
        {
            ReadCigarData next = readCigars.get(i);

            if(next.Elements.size() != first.Elements.size())
            {
                hasMismatch = true;
                break;
            }

            for(int j = 1; j < first.Elements.size(); ++j) // ignore first M
            {
                if(first.Elements.get(j).Operator != next.Elements.get(j).Operator)
                {
                    hasMismatch = true;
                    break;
                }

                if(first.Elements.get(j).PosStart != next.Elements.get(j).PosStart)
                {
                    hasMismatch = true;
                    break;
                }
            }

            if(hasMismatch)
                break;
        }

        if(!hasMismatch)
        {
            mBaseBuilder.buildReadBases(reads, consensusState);
            consensusState.setOutcome(INDEL_MATCH);
            buildIndelCigar(consensusState,reads.get(0));
        }

        /*
        // check for differences and incompatibilities
        int firstElementCount = readCigars.get(0).size();

        for(int i = 0; i < reads.size(); ++i)
        {
            if(readCigars.get(i).size() != firstElementCount)
            {
                consensusState.setOutcome(INDEL_MISMATCH);
                return;
            }
        }

        List<CigarElementData> firstElementData = readCigars.get(0);

        for(int i = 0; i < firstElementData.size(); ++i)
        {
            CigarElementData firstElement = firstElementData.get(i);

            for(int j = 1; j < reads.size(); ++j)
            {
                CigarElementData otherElement = readCigars.get(j).get(i);

                if(otherElement.Operator != firstElement.Operator)
                {
                    consensusState.setOutcome(INDEL_MISMATCH);
                    return;
                }
            }
        }
        */

        // List<CigarElementData> consensusCigarElements = Lists.newArrayListWithExpectedSize(reads.size());


    }

    private class ReadCigarData
    {
        public boolean Compatible;

        public final List<CigarElementData> Elements;

        public ReadCigarData(final List<CigarElementData> elements)
        {
            Elements = elements;
            Compatible = true;
        }
    }

    private class CigarElementData
    {
        public final int ReadIndexStart;
        public final int PosStart;
        public final CigarOperator Operator;
        public final int Length;

        public CigarElementData(final int readIndexStart, final int posStart, final CigarOperator operator, final int length)
        {
            ReadIndexStart = readIndexStart;
            PosStart = posStart;
            Operator = operator;
            Length = length;
        }

        public String toString() { return format("i(%d) pos(%d) e(%d%s)", ReadIndexStart, PosStart, Length, Operator); }
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

    private static final byte INVALID_BASE = -1;

    private class ReadState
    {
        public final SAMRecord Read;

        private int mReadIndex;
        private int mCigarIndex;
        private int mCigarElementIndex;
        private int mCigarElementLength;
        private boolean mExhausted;

        public ReadState(final SAMRecord read)
        {
            Read = read;

            if(!read.getReadNegativeStrandFlag())
            {
                mReadIndex = 0;
                mCigarIndex = 0;
            }
            else
            {
                mReadIndex = read.getReadBases().length - 1;
                mCigarIndex = read.getCigar().getCigarElements().size() - 1;
            }

            mCigarElementLength = read.getCigar().getCigarElement(mCigarIndex).getLength();
            mCigarElementIndex = 0;
            mExhausted = false;
        }

        public boolean exhausted() { return mExhausted; }

        public void moveNext(boolean isForward)
        {
            if(mExhausted)
                return;

            if(isForward)
            {
                ++mReadIndex;

                if(mReadIndex >= Read.getReadBases().length)
                {
                    mExhausted = true;
                    return;
                }
            }
            else
            {
                --mReadIndex;

                if(mReadIndex < 0)
                {
                    mExhausted = true;
                    return;
                }
            }

            ++mCigarElementIndex;

            if(mCigarElementIndex >= mCigarElementLength)
            {
                if(isForward)
                    ++mCigarIndex;
                else
                    --mCigarIndex;

                mCigarElementLength = Read.getCigar().getCigarElement(mCigarIndex).getLength();
                mCigarElementIndex = 0;
            }
        }

        public byte currentBase() { return mExhausted ? INVALID_BASE : Read.getReadBases()[mReadIndex]; }
        public byte currentBaseQual() { return mExhausted ? INVALID_BASE : Read.getBaseQualities()[mReadIndex]; }
        public CigarOperator currentCigarType() { return Read.getCigar().getCigarElement(mCigarIndex).getOperator(); }

        public String toString()
        {
            return format("index(%d) cigar(%d element=%d/%d) %s",
                    mReadIndex, mCigarIndex, mCigarElementIndex, mCigarElementLength, mExhausted ? "exhausted" : "active");
        }
    }

}
