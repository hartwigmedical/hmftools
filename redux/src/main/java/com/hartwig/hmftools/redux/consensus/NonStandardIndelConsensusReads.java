package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_ERROR_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.DUPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.SIMPLEX_QUAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;

import static htsjdk.samtools.CigarOperator.M;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.apache.commons.lang3.NotImplementedException;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public abstract class NonStandardIndelConsensusReads
{
    private final static byte ANY_BASE = (byte) 'N';

    protected final RefGenome mRefGenome;
    protected int mChromosomeLength;

    public NonStandardIndelConsensusReads(final RefGenome refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLength = 0;
    }

    public abstract void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState);

    public void setChromosomeLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }

    public static NonStandardIndelConsensusReads fromSequencingType(final SequencingType sequencingType, final RefGenome refGenome)
    {
        if(sequencingType == SBX)
            return new SBXBaseBuilder(refGenome);

        return null;
    }

    private static List<CigarElement> normalizeCigar(final List<CigarElement> cigar)
    {
        List<CigarElement> normalizedElements = Lists.newArrayList();
        for(CigarElement el : cigar)
        {
            switch(el.getOperator())
            {
                case H:
                    break;
                case S:
                    normalizedElements.add(el);
                    break;
                case M:
                case EQ:
                case X:
                    normalizedElements.add(new CigarElement(el.getLength(), M));
                    break;
                default:
                    // TODO:
                    throw new NotImplementedException("TODO");
            }
        }

        return normalizedElements;
    }

    private static List<CigarElement> mergeRepeatedCigarOps(final List<CigarElement> elements)
    {
        List<CigarElement> mergedElements = Lists.newArrayList();
        CigarElement lastElement = null;
        for(CigarElement el : elements)
        {
            if(lastElement == null)
            {
                mergedElements.add(el);
                lastElement = el;
                continue;
            }

            if(lastElement.getOperator() == el.getOperator())
            {
                lastElement = new CigarElement(lastElement.getLength() + el.getLength(), lastElement.getOperator());
                mergedElements.set(mergedElements.size() - 1, lastElement);
                continue;
            }

            mergedElements.add(el);
            lastElement = el;
        }

        return mergedElements;
    }

    private static class TrimmedCigar
    {
        private final String mReadName;
        private final int mUnclippedFivePrimePos;
        private final List<CigarElement> mForwardElements;
        private final String mForwardCigarStr;
        private final int mAlignmentLength;

        public TrimmedCigar(final String readName, int unclippedFivePrimePos, final List<CigarElement> forwardElements)
        {
            mReadName = readName;
            mUnclippedFivePrimePos = unclippedFivePrimePos;
            mForwardElements = forwardElements;
            mForwardCigarStr = forwardElements.stream().map(CigarElement::toString).collect(Collectors.joining());
            mAlignmentLength = forwardElements.stream().mapToInt(x -> x.getOperator() == M ? x.getLength() : 0).sum();
        }

        public String readName()
        {
            return mReadName;
        }

        public int unclippedFivePrimePos()
        {
            return mUnclippedFivePrimePos;
        }

        public List<CigarElement> forwardElements()
        {
            return mForwardElements;
        }

        public String forwardCigarStr()
        {
            return mForwardCigarStr;
        }

        public int alignmentLength()
        {
            return mAlignmentLength;
        }
    }

    private static class CigarInfo
    {
        private final String mReadName;
        private final int mUnclippedFivePrimePos;
        private final Deque<CigarElement> mForwardElements;

        private int mRemainingReadLength;

        public CigarInfo(final String readName, int unclippedFivePrimePos, final Collection<CigarElement> forwardElements)
        {
            mReadName = readName;
            mUnclippedFivePrimePos = unclippedFivePrimePos;
            mForwardElements = new ArrayDeque<>(forwardElements);
            mRemainingReadLength = forwardElements.stream().mapToInt(CigarElement::getLength).sum();
        }

        public TrimmedCigar trim(int readLength)
        {
            List<CigarElement> head = Lists.newArrayList();
            int lengthRemaining = readLength;
            while(lengthRemaining > 0)
            {
                CigarElement el = mForwardElements.pop();
                if(el.getLength() <= lengthRemaining)
                {
                    head.add(el);
                    lengthRemaining -= el.getLength();
                    continue;
                }

                int elLengthRemaining = el.getLength() - lengthRemaining;
                CigarOperator op = el.getOperator();
                mForwardElements.push(new CigarElement(elLengthRemaining, op));
                head.add(new CigarElement(lengthRemaining, op));
                lengthRemaining = 0;
            }

            return new TrimmedCigar(mReadName, mUnclippedFivePrimePos, head);
        }

        public int remainingReadLength()
        {
            return mRemainingReadLength;
        }

        public boolean isExhausted()
        {
            return mRemainingReadLength == 0;
        }
    }

    private static TrimmedCigar getTrimmedConsensusCigar(final List<TrimmedCigar> trimmedCigars)
    {
        Map<String, Integer> cigarFreqs = Maps.newHashMap();
        int maxFreq = 0;
        for(TrimmedCigar trimmedCigar : trimmedCigars)
        {
            String cigarStr = trimmedCigar.forwardCigarStr();
            cigarFreqs.merge(cigarStr, 1, Integer::sum);
            maxFreq = max(maxFreq, cigarFreqs.get(cigarStr));
        }

        Set<String> maxFreqCigarStrs = Sets.newHashSet();
        for(Map.Entry<String, Integer> cigarStrAndFreq : cigarFreqs.entrySet())
        {
            String cigarStr = cigarStrAndFreq.getKey();
            int freq = cigarStrAndFreq.getValue();
            if(freq == maxFreq)
                maxFreqCigarStrs.add(cigarStr);
        }

        List<TrimmedCigar> maxFreqTrimmedCigars = trimmedCigars.stream().filter(x -> maxFreqCigarStrs.contains(x.forwardCigarStr())).collect(Collectors.toList());
        return Collections.min(maxFreqTrimmedCigars, (final TrimmedCigar o1, final TrimmedCigar o2) -> {
            if(o1.alignmentLength() != o2.alignmentLength())
                return o2.alignmentLength() - o1.alignmentLength();

            return o1.readName().compareTo(o2.readName());
        });
    }

    private static class ConsensusCigarInfo
    {
        public final int UnclippedFivePrimePos;
        public final List<CigarElement> CigarElements;

        public ConsensusCigarInfo(final int unclippedFivePrimePos, final List<CigarElement> cigarElements)
        {
            UnclippedFivePrimePos = unclippedFivePrimePos;
            CigarElements = cigarElements;
        }
    }

    private static ConsensusCigarInfo buildConsensusCigar(final List<SAMRecord> reads, boolean isForward)
    {
        List<CigarInfo> cigarInfos = Lists.newArrayList();
        for(SAMRecord read : reads)
        {
            int unclippedFivePrimePos = read.getReadNegativeStrandFlag() ? read.getUnclippedEnd() : read.getUnclippedStart();
            List<CigarElement> cigar = normalizeCigar(read.getCigar().getCigarElements());
            if(!isForward)
                Collections.reverse(cigar);

            cigarInfos.add(new CigarInfo(read.getReadName(), unclippedFivePrimePos, cigar));
        }

        List<TrimmedCigar> consensus = Lists.newArrayList();

        // TODO: while(true)
        while(true)
        {
            int minReadLength = cigarInfos.stream().filter(x -> !x.isExhausted()).mapToInt(x -> x.remainingReadLength()).min().orElse(0);
            if(minReadLength == 0)
                break;

            List<TrimmedCigar> trimmedCigars = cigarInfos.stream().filter(x -> !x.isExhausted()).map(x -> x.trim(minReadLength)).collect(Collectors.toList());
            consensus.add(getTrimmedConsensusCigar(trimmedCigars));
        }

        int unclippedFivePrimePos = consensus.get(0).unclippedFivePrimePos();
        List<CigarElement> consensusCigar = Lists.newArrayList();
        for(TrimmedCigar trimmedCigar : consensus)
            consensusCigar.addAll(trimmedCigar.forwardElements());

        consensusCigar = mergeRepeatedCigarOps(consensusCigar);

        if(!isForward)
            Collections.reverse(consensusCigar);

        return new ConsensusCigarInfo(unclippedFivePrimePos, consensusCigar);
    }

    private static void updateConsensusStateBoundaries(int unclippedFivePrimePos, final ConsensusState consensusState)
    {
        int unclippedStart;
        int unclippedEnd;
        if(consensusState.IsForward)
        {
            unclippedStart = unclippedFivePrimePos;
            // TODO: move into util function
            unclippedEnd = unclippedStart + consensusState.CigarElements.stream().mapToInt(x -> x.getLength()).sum();
        }
        else
        {
            unclippedEnd = unclippedFivePrimePos;
            unclippedStart = unclippedEnd - consensusState.CigarElements.stream().mapToInt(x -> x.getLength()).sum();
        }

        int readStart = unclippedStart + leftSoftClipLength(consensusState.CigarElements);
        int readEnd = unclippedEnd - rightSoftClipLength(consensusState.CigarElements);
        consensusState.setBoundaries(unclippedStart, unclippedEnd, readStart, readEnd);
    }

    public static class SBXBaseBuilder extends NonStandardIndelConsensusReads
    {
        public SBXBaseBuilder(final RefGenome refGenome)
        {
            super(refGenome);
        }

        @Override
        public void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState)
        {
            ConsensusCigarInfo consensusCigarInfo = buildConsensusCigar(reads, consensusState.IsForward);
            int consensusReadLength = consensusCigarInfo.CigarElements.stream().mapToInt(CigarElement::getLength).sum();
            consensusState.setBaseLength(consensusReadLength);
            buildReadBases(reads, consensusCigarInfo.UnclippedFivePrimePos, consensusState);
            consensusState.CigarElements.addAll(consensusCigarInfo.CigarElements);
            updateConsensusStateBoundaries(consensusCigarInfo.UnclippedFivePrimePos, consensusState);
        }

        private void buildReadBases(final List<SAMRecord> reads, int unclippedFivePrimePos, final ConsensusState consensusState)
        {
            String chromosome = reads.get(0).getReferenceName();
            mChromosomeLength = mRefGenome.getChromosomeLength(chromosome);

            int baseLength = consensusState.Bases.length;
            int readCount = reads.size();

            byte[] locationBases = new byte[readCount];
            byte[] locationQuals = new byte[readCount];

            for(int i = 0; i < baseLength; ++i)
            {
                for(int r = 0; r < readCount; ++r)
                {
                    SAMRecord read = reads.get(r);
                    int baseIndex = consensusState.IsForward ? i : read.getReadLength() - 1 - i;
                    if(baseIndex < 0 || baseIndex >= read.getReadLength())
                    {
                        locationBases[r] = NO_BASE;
                        continue;
                    }

                    locationBases[r] = read.getReadBases()[baseIndex];
                    locationQuals[r] = read.getBaseQualities()[baseIndex];
                }

                int basePosition = consensusState.IsForward ? unclippedFivePrimePos + i : unclippedFivePrimePos - i;

                if(basePosition < 1 || basePosition > mChromosomeLength)
                    basePosition = INVALID_POSITION; // protect against over-runs from soft-clips - rare but possible

                byte[] consensusBaseAndQual = determineBaseAndQual(locationBases, locationQuals, chromosome, basePosition);
                int baseIndex = consensusState.IsForward ? i : baseLength - 1 - i;
                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }
        }

        private byte[] determineBaseAndQual(final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position)
        {
            int baseCount = locationBases.length;
            if(baseCount == 1)
            {
                boolean isQualZero = locationQuals[0] == (byte) 0;
                if(isQualZero && position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new byte[] { refBase, (byte) 1 };
                }

                return new byte[] { locationBases[0], locationQuals[0] };
            }

            Map<Byte, int[]> baseCountsByQual = Maps.newHashMap();
            for(int i = 0; i < baseCount; i++)
            {
                if(locationBases[i] == NO_BASE)
                    continue;

                byte base = locationBases[i];
                int baseIdx = baseIndex(base);
                if(baseIdx < 0 || baseIdx >= DNA_BASE_BYTES.length)
                    continue;

                byte qual = locationQuals[i];
                baseCountsByQual.computeIfAbsent(qual, key -> new int[] { 0, 0, 0, 0 });
                int[] baseCounts = baseCountsByQual.get(qual);
                baseCounts[baseIdx]++;
            }

            if(baseCountsByQual.isEmpty())
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new byte[] { refBase, (byte) 1 };
                }

                return new byte[] { ANY_BASE, (byte) 0 };
            }

            if(!baseCountsByQual.containsKey((byte) SIMPLEX_QUAL) && !baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new byte[] { refBase, (byte) 1 };
                }

                return new byte[] { ANY_BASE, (byte) 0 };
            }

            if(!baseCountsByQual.containsKey((byte) DUPLEX_QUAL))
            {
                int[] simplexCounts = baseCountsByQual.get((byte) SIMPLEX_QUAL);
                int maxIdx = -1;
                int maxCount = -1;
                int totalCount = 0;
                for(int i = 0; i < simplexCounts.length; i++)
                {
                    totalCount += simplexCounts[i];
                    if(simplexCounts[i] > maxCount)
                    {
                        maxCount = simplexCounts[i];
                        maxIdx = i;
                    }
                }

                if(2 * maxCount > totalCount)
                {
                    byte base = DNA_BASE_BYTES[maxIdx];
                    byte qual = (byte) SIMPLEX_QUAL;
                    return new byte[] { base, qual };
                }

                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new byte[] { refBase, (byte) 1 };
                }

                return new byte[] { ANY_BASE, (byte) 0 };
            }

            int[] duplexCounts = baseCountsByQual.get((byte) DUPLEX_QUAL);
            int maxIdx = -1;
            boolean multipleMax = false;
            int maxCount = -1;
            int totalCount = 0;
            for(int i = 0; i < duplexCounts.length; i++)
            {
                totalCount += duplexCounts[i];
                if(duplexCounts[i] > maxCount)
                {
                    maxCount = duplexCounts[i];
                    maxIdx = i;
                    multipleMax = false;
                }
                else if(duplexCounts[i] == maxCount)
                    multipleMax = true;
            }

            int[] duplexErrorCounts = baseCountsByQual.get((byte) DUPLEX_ERROR_QUAL);
            if(duplexErrorCounts != null)
            {
                for(int duplexErrorCount : duplexErrorCounts)
                    totalCount += duplexErrorCount;
            }

            if(2 * maxCount > totalCount)
            {
                byte base = DNA_BASE_BYTES[maxIdx];
                byte qual = (byte) DUPLEX_QUAL;
                return new byte[] { base, qual };
            }

            if(multipleMax)
            {
                if(position != INVALID_POSITION)
                {
                    byte refBase = mRefGenome.getRefBase(chromosome, position);
                    return new byte[] { refBase, (byte) 1 };
                }

                return new byte[] { ANY_BASE, (byte) 0 };
            }

            byte base = DNA_BASE_BYTES[maxIdx];
            byte qual = (byte) SIMPLEX_QUAL;
            return new byte[] { base, qual };
        }
    }
}
