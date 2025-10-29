package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_REPEAT_BASE_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.DNA_BASE_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffInfo.UNSET;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.MATCH;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.hasMinRepeats;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHighBaseQual;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public class SequenceBuilder
{
    private final boolean mBuildForwards;
    private final int mInitialBaseLength;
    private final List<ReadParseState> mReads;

    private int mCurrentIndex;
    private byte[] mBases;
    private byte[] mBaseQuals;

    private final List<RepeatInfo> mRepeats;

    public SequenceBuilder(final List<ReadParseState> reads, boolean buildForwards, int initialBaseLength)
    {
        mBuildForwards = buildForwards;
        mInitialBaseLength = initialBaseLength;
        mReads = reads;

        mBases = new byte[initialBaseLength];
        mBaseQuals = new byte[initialBaseLength];
        mCurrentIndex = -1;
        mRepeats = Lists.newArrayList();

        buildSequence();
        trimFinalSequence();
    }

    public byte[] bases() { return mBases; }
    public String baseString() { return new String(mBases); }
    public byte[] baseQuals() { return mBaseQuals; }
    public List<ReadParseState> reads() { return mReads; }
    public List<RepeatInfo> repeats() { return mRepeats; }

    private void buildSequence()
    {
        mCurrentIndex = mBuildForwards ? 0 : mBases.length - 1;

        List<ReadParseState> activeReads = Lists.newArrayList(mReads);

        while(mCurrentIndex >= 0 && mCurrentIndex < mBases.length)
        {
            byte consensusBase = 0;
            byte consensusMaxQual = 0;
            int consensusHighQualCount = 0;
            int consensusMedQualCount = 0;

            byte[] maxQuals = null; // max qual per base
            int[] highQuals = null; // count of reads with high-quals per base
            int[] mediumQuals = null; // count of reads with high-quals per base

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);

                /*
                if(read.exceedsMaxMismatches()) // relevant when building the final consensus
                    continue;
                */

                byte base = read.currentBase();
                byte qual = read.currentQual();
                boolean isHighQual = isHighBaseQual(qual);

                if(Nucleotides.baseIndex(base) < 0)
                {
                    base = DNA_N_BYTE;
                    qual = 0;
                    isHighQual = false;
                }

                boolean hasMismatch = maxQuals != null;

                if(maxQuals == null)
                {
                    if(consensusBase == 0 || (base != consensusBase && belowMinQual(consensusMaxQual) && aboveMinQual(qual)))
                    {
                        // set first or replace with first high qual
                        consensusBase = base;
                        consensusMaxQual = qual;

                        if(isHighQual)
                            ++consensusHighQualCount;
                        else
                            ++consensusMedQualCount;
                    }
                    else if(base == consensusBase)
                    {
                        consensusMaxQual = maxQual(qual, consensusMaxQual);

                        if(isHighQual)
                            ++consensusHighQualCount;
                        else
                            ++consensusMedQualCount;
                    }
                    else if(base != consensusBase && belowMinQual(qual))
                    {
                        // low-qual disagreement - ignore regardless of consensus qual
                    }
                    else
                    {
                        hasMismatch = true;

                        // high-qual mismatch so start tracking frequencies for each base
                        maxQuals = new byte[DNA_BASE_COUNT];
                        highQuals = new int[DNA_BASE_COUNT];
                        mediumQuals = new int[DNA_BASE_COUNT];

                        // back port existing counts to the per-base arrays
                        int baseIndex = Nucleotides.baseIndex(consensusBase);
                        maxQuals[baseIndex] = consensusMaxQual;
                        highQuals[baseIndex] = consensusHighQualCount;
                        mediumQuals[baseIndex] = consensusMedQualCount;
                    }
                }

                if(hasMismatch)
                {
                    int baseIndex = Nucleotides.baseIndex(base);

                    maxQuals[baseIndex] = maxQual(maxQuals[baseIndex], qual);

                    if(isHighQual)
                        ++highQuals[baseIndex];
                    else
                        ++mediumQuals[baseIndex];
                }
            }

            if(maxQuals != null)
            {
                // take the bases with high number of high-qual bases, otherwise higher number of mediums
                int maxHighQual = 0;
                int maxMediumQual = 0;
                int maxBaseIndex = 0;

                for(int b = 0; b < maxQuals.length; ++b)
                {
                    if(highQuals[b] > maxHighQual)
                    {
                        maxHighQual = highQuals[b];
                        maxMediumQual = mediumQuals[b];
                        maxBaseIndex = b;
                    }
                    else if(highQuals[b] == maxHighQual && mediumQuals[b] > maxMediumQual)
                    {
                        maxHighQual = highQuals[b];
                        maxMediumQual = mediumQuals[b];
                        maxBaseIndex = b;
                    }
                }

                consensusBase = maxBaseIndex < DNA_BASE_BYTES.length ? DNA_BASE_BYTES[maxBaseIndex] : DNA_N_BYTE;
                consensusMaxQual = maxQuals[maxBaseIndex];

                SequenceDiffInfo[] seqDiffInfos = assessDifferences(consensusBase, consensusMaxQual, activeReads);

                applyReadMismatches(activeReads, consensusBase, consensusMaxQual, seqDiffInfos);
            }
            else
            {
                mBases[mCurrentIndex] = consensusBase;
                mBaseQuals[mCurrentIndex] = consensusMaxQual;

                mCurrentIndex += mBuildForwards ? 1 : -1;

                activeReads.forEach(x -> x.moveNext());
            }

            // purge exhausted reads
            int readIndex = 0;
            while(readIndex < activeReads.size())
            {
                ReadParseState read = activeReads.get(readIndex);

                if(read.exhausted() || !read.isValid())
                    activeReads.remove(readIndex);
                else
                    ++readIndex;
            }

            if(activeReads.isEmpty())
                break;
        }
    }

    private void applyReadMismatches(
            final List<ReadParseState> activeReads, byte consensusBase, byte consensusQual, final SequenceDiffInfo[] seqDiffInfos)
    {
        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);

            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];

            read.moveOnMatchType(seqDiffInfo);

            if(seqDiffInfo.Type != MATCH)
                read.addMismatchInfo(seqDiffInfo);
        }

        RepeatInfo repeatInfo = Arrays.stream(seqDiffInfos).anyMatch(x -> x.RepeatCount > 0) ? mRepeats.get(mRepeats.size() - 1) : null;

        if(repeatInfo == null)
        {
            mBases[mCurrentIndex] = consensusBase;
            mBaseQuals[mCurrentIndex] = consensusQual;
            mCurrentIndex += mBuildForwards ? 1 : -1;
            return;
        }

        int repeatLength = repeatInfo.baseLength();

        // adjust each read's index relative to the consensus repeat
        int previousReadBases = abs(mCurrentIndex - repeatInfo.Index);
        int remainingRepeatBases = repeatInfo.length() - previousReadBases;

        int repeatBaseIndex = mBuildForwards ? 0 : repeatLength - 1;
        String repeatBases = repeatInfo.Bases;

        for(int i = 0; i < remainingRepeatBases; ++i)
        {
            // apply the repeat to the consensus bases
            byte repeatBase = (byte)repeatBases.charAt(repeatBaseIndex);
            mBases[mCurrentIndex] = repeatBase;
            mBaseQuals[mCurrentIndex] = consensusQual;

            if(mBuildForwards)
            {
                ++mCurrentIndex;
                ++repeatBaseIndex;

                if(repeatBaseIndex == repeatLength)
                    repeatBaseIndex = 0;
            }
            else
            {
                --mCurrentIndex;
                --repeatBaseIndex;

                if(repeatBaseIndex < 0)
                    repeatBaseIndex = repeatLength - 1;
            }
        }

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);

            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];
            int remainingReadRepeatBases = seqDiffInfo.RepeatCount * repeatInfo.baseLength() - previousReadBases;
            read.moveNextBases(remainingReadRepeatBases);
        }
    }

    private static final int NEXT_BASE_CHECK_COUNT = 3;

    private SequenceDiffInfo[] assessDifferences(byte consensusBase, byte consensusQual, final List<ReadParseState> activeReads)
    {
        // assess each read relative to the consensus base to establish if an indel or repeat difference is the cause
        SequenceDiffInfo[] seqDiffInfos = new SequenceDiffInfo[activeReads.size()];

        for(int i = 0; i < seqDiffInfos.length; ++i)
        {
            seqDiffInfos[i] = UNSET;
        }

        int readCount = activeReads.size();

        // first check if bases match after this single base mismatch, ie no indel-based differences
        List<byte[]> allReadNextBases = Lists.newArrayListWithCapacity(readCount);
        List<byte[]> consensusReadNextBases = Lists.newArrayListWithCapacity(readCount);

        SequenceDiffInfo matchInfo = SequenceDiffInfo.fromMatch(mCurrentIndex);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);

            byte[] readNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);

            allReadNextBases.add(readNextBases); // stored even if null

            boolean baseMatch = basesMatch(consensusBase, read.currentBase(), consensusQual, read.currentQual());

            if(read.currentBase() == consensusBase)
                consensusReadNextBases.add(readNextBases);

            if(baseMatch)
            {
                seqDiffInfos[readIndex] = matchInfo;
            }
            else if(readNextBases == null)
            {
                seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(mCurrentIndex, read.currentBase());
            }
        }

        // first check if the consensus base matches the prior base for the reads which match
        byte[] consensusNextBases = findConsensusNextBases(consensusReadNextBases);

        if(consensusNextBases != null)
        {
            int postConsensusMismatchCount = 0;

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                if(seqDiffInfos[readIndex].Type == MATCH)
                    continue;

                ReadParseState read = activeReads.get(readIndex);
                byte[] readNextBases = allReadNextBases.get(readIndex);

                if(readNextBases == null)
                    continue;

                if(basesMatch(readNextBases, consensusNextBases))
                {
                    seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(mCurrentIndex, read.currentBase());
                }
                else
                {
                    ++postConsensusMismatchCount;
                }
            }

            if(postConsensusMismatchCount == 0)
            {
                // all mismatches are explained by a SNV but subsequent bases match
                return seqDiffInfos;
            }

            assessNovelIndelDifferences(consensusBase, consensusNextBases, activeReads, allReadNextBases, seqDiffInfos);

            int remainingMismatchCount = (int)Arrays.stream(seqDiffInfos).filter(x -> x == UNSET).count();

            if(remainingMismatchCount == 0)
                return seqDiffInfos;
        }

        // otherwise check for a repeat expansion / contraction
        assessRepeatDifferences(activeReads, seqDiffInfos);

        return seqDiffInfos;
    }

    private void assessNovelIndelDifferences(
            byte consensusBase, final byte[] consensusNextBases, final List<ReadParseState> activeReads,
            final List<byte[]> allReadNextBases, final SequenceDiffInfo[] seqDiffInfos)
    {
        // currently only searches for 1-base indels which aren't repeat expansions or contractions
        // non-consensus read (NCR) current + next X-1 bases match CR subsequent bases = 1-base delete
        // NCR next bases match CR current + subsequent X-1 bases = 1-base ins
        byte[] consensusCurrentAndNext = currentAndNextBase(consensusBase, consensusNextBases);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);
            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];
            byte[] readNextBases = allReadNextBases.get(readIndex);

            if(seqDiffInfo.Type == MATCH || seqDiffInfo.Type == BASE || readNextBases == null)
                continue;

            if(basesMatch(consensusCurrentAndNext, readNextBases))
            {
                // read has inserted the current base
                seqDiffInfos[readIndex] = new SequenceDiffInfo(
                        mCurrentIndex, String.valueOf((char)read.currentBase()), SequenceDiffType.INSERT, 0);
            }
            else
            {
                byte[] readCurrentAndNext = currentAndNextBase(read.currentBase(), readNextBases);
                if(basesMatch(consensusNextBases, readCurrentAndNext))
                {
                    seqDiffInfos[readIndex] = new SequenceDiffInfo(mCurrentIndex, "", SequenceDiffType.DELETE, 0);
                }
            }
        }
    }

    private void assessRepeatDifferences(final List<ReadParseState> activeReads, final SequenceDiffInfo[] seqDiffInfos)
    {
        // assume that the repeat ends for at least some of the reads at the base prior to the mismatch base,
        // so search for repeats around this location

        // first of all look for single-base expansions since this is the most common
        String maxRepeatBases = null;
        int maxReadsWithRepeatCount = 0; // number of reads with evidence of the repeat
        int[] maxRepeatReadRepeatCounts = null; // the repeat counts per read for the most likely repeat
        int maxRepeatPreviousRepeatCount = 0;

        for(int repeatLength = 1; repeatLength <= MAX_REPEAT_BASE_COUNT; ++repeatLength)
        {
            String repeatBases = null;

            if(repeatLength == 1)
            {
                byte previousBase = mBuildForwards ? mBases[mCurrentIndex - 1] : mBases[mCurrentIndex + 1];
                repeatBases = String.valueOf((char)previousBase);
            }
            else
            {
                if(mBuildForwards)
                {
                    if(mCurrentIndex - repeatLength < 0)
                        break;

                    repeatBases = new String(Arrays.copyOfRange(mBases, mCurrentIndex - repeatLength, mCurrentIndex));
                }
                else
                {
                    if(mCurrentIndex + repeatLength >= mBases.length)
                        break;

                    repeatBases = new String(Arrays.copyOfRange(mBases, mCurrentIndex + 1, mCurrentIndex + repeatLength + 1));
                }
            }

            int repeatIndexStart = mBuildForwards ? mCurrentIndex - 1 : mCurrentIndex + 1;
            int previousRepeatCount = RepeatInfo.getRepeatCount(mBases, repeatBases, repeatIndexStart, !mBuildForwards);

            if(!hasMinRepeats(repeatLength, previousRepeatCount))
                continue;

            int[] repeatReadRepeatCounts = new int[activeReads.size()];
            int nonZeroRepeatReadCount = 0;

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);

                int repeatCount = getRepeatCount(read, repeatBases, previousRepeatCount);
                repeatReadRepeatCounts[readIndex] = repeatCount;

                if(repeatCount > 0)
                    ++nonZeroRepeatReadCount;
            }

            if(nonZeroRepeatReadCount > maxReadsWithRepeatCount)
            {
                maxReadsWithRepeatCount = nonZeroRepeatReadCount;
                maxRepeatReadRepeatCounts = repeatReadRepeatCounts;
                maxRepeatBases = repeatBases;
                maxRepeatPreviousRepeatCount = previousRepeatCount;

                if(nonZeroRepeatReadCount == activeReads.size())
                    break;
            }
        }

        if(maxRepeatBases == null)
            return;

        // repeat count must meet minimums to be registered
        int maxObservedRepeatCount = Arrays.stream(maxRepeatReadRepeatCounts).max().orElse(0);

        if(!hasMinRepeats(maxRepeatBases.length(), maxObservedRepeatCount))
            return;

        // must impact more than half the reads
        if(maxReadsWithRepeatCount <= 0.5 * activeReads.size())
            return;

        int consensusRepeatCount = findConsensusRepeatCount(seqDiffInfos, maxRepeatReadRepeatCounts);

        // register this repeat against the consensus sequence
        int consensusRepeatOffset = maxRepeatPreviousRepeatCount * maxRepeatBases.length();
        int consensusRepeatIndexStart = mBuildForwards ? mCurrentIndex - consensusRepeatOffset : mCurrentIndex + consensusRepeatOffset;

        RepeatInfo repeatInfo = new RepeatInfo(consensusRepeatIndexStart, maxRepeatBases, consensusRepeatCount);
        mRepeats.add(repeatInfo);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];
            int readRepeatCount = maxRepeatReadRepeatCounts[readIndex];

            SequenceDiffType diffType = readRepeatCount == consensusRepeatCount ? MATCH : SequenceDiffType.REPEAT;

            seqDiffInfos[readIndex] = new SequenceDiffInfo(mCurrentIndex, maxRepeatBases, diffType, readRepeatCount);
        }
    }

    private int getRepeatCount(final ReadParseState read, final String repeat, int previousRepeatCount)
    {
        // search forward from the current index
        int countForward = max(RepeatInfo.getRepeatCount(read.read().getBases(), repeat, read.readIndex(), mBuildForwards), 0);
        return previousRepeatCount + countForward;
    }

    private static int findConsensusRepeatCount(final SequenceDiffInfo[] seqDiffInfos, final int[] readRepeatCounts)
    {
        Map<Integer,Integer> repeatFrequencies = Maps.newHashMap();

        for(int readIndex = 0; readIndex < readRepeatCounts.length; ++readIndex)
        {
            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];

            if(seqDiffInfo.Type == MATCH)
            {
                int readRepeatCount = readRepeatCounts[readIndex];

                Integer readCount = repeatFrequencies.get(readRepeatCount);
                repeatFrequencies.put(readRepeatCount, readCount == null ? 1 : readCount + 1);
            }
        }

        int maxRepeatFreq = 0;
        int maxRepeatCount = 0;

        for(Map.Entry<Integer,Integer> entry : repeatFrequencies.entrySet())
        {
            if(entry.getValue() > maxRepeatFreq)
            {
                maxRepeatCount = entry.getKey();
                maxRepeatFreq = entry.getValue();
            }
        }

        return maxRepeatCount;
    }

    private byte[] currentAndNextBase(byte currentBase, final byte[] nextBases)
    {
        byte[] currentAndNext = new byte[nextBases.length];

        if(mBuildForwards)
        {
            int index = 0;
            currentAndNext[index++] = currentBase;

            for(int i = 0; i < nextBases.length - 1; ++i)
            {
                currentAndNext[index++] = nextBases[i];
            }
        }
        else
        {
            int index = nextBases.length - 1;
            currentAndNext[index--] = currentBase;

            for(int i = nextBases.length - 1; i >= 1; --i)
            {
                currentAndNext[index--] = nextBases[i];
            }
        }

        return currentAndNext;
    }

    private static byte[] findConsensusNextBases(final List<byte[]> consensusReadNextBases)
    {
        // find the following set of consensus bases
        byte[] consensusNextBases = new byte[NEXT_BASE_CHECK_COUNT];

        Map<byte[],Integer> consensusNextBasesFreq = null;

        int consensusReadNextBaseMatchCount = 0;
        int consensusReadCount = consensusReadNextBases.size();

        for(byte[] readNextBases : consensusReadNextBases)
        {
            if(readNextBases == null) // cannot be assessed
                continue;

            if(consensusNextBases == null)
            {
                consensusNextBases = readNextBases;
                consensusReadNextBaseMatchCount = 1;
            }
            else if(consensusNextBasesFreq == null)
            {
                if(basesMatch(readNextBases, consensusNextBases))
                {
                    ++consensusReadNextBaseMatchCount;
                }
                else
                {
                    consensusNextBasesFreq = Maps.newHashMap();
                    consensusNextBasesFreq.put(consensusNextBases, consensusReadNextBaseMatchCount);

                    consensusNextBasesFreq.put(readNextBases, 1);
                }
            }
            else
            {
                boolean found = false;

                for(Map.Entry<byte[],Integer> entry : consensusNextBasesFreq.entrySet())
                {
                    if(basesMatch(readNextBases, entry.getKey()))
                    {
                        entry.setValue(entry.getValue() + 1);
                        found = true;
                        break;
                    }
                }

                if(!found)
                    consensusNextBasesFreq.put(readNextBases, 1);
            }
        }

        if(consensusNextBasesFreq != null)
        {
            for(Map.Entry<byte[], Integer> entry : consensusNextBasesFreq.entrySet())
            {
                if(entry.getValue() > consensusReadNextBaseMatchCount)
                {
                    consensusReadNextBaseMatchCount = entry.getValue();
                    consensusNextBases = entry.getKey();
                }
            }
        }

        // only take the next consensus bases if they are present in the majority of consensus reads
        double matchedFraction = consensusReadNextBaseMatchCount / (double)consensusReadCount;

        return matchedFraction > 0.5 ? consensusNextBases : null;
    }

    private void trimFinalSequence()
    {
        int currentIndex = mBuildForwards ? 0 : mBases.length - 1;
        int validLength = 0;

        while(currentIndex >= 0 && currentIndex < mBases.length)
        {
            if(mBases[currentIndex] == 0)
                break;

            ++validLength;
            currentIndex += mBuildForwards ? 1 : -1;
        }

        if(validLength == mBases.length)
            return;

        int reduction = mBases.length - validLength;
        int startIndex = mBuildForwards ? 0 : reduction;
        int endIndex = startIndex + validLength - 1;
        mBases = subsetArray(mBases, startIndex, endIndex);
        mBaseQuals = subsetArray(mBaseQuals, startIndex, endIndex);
    }
}
