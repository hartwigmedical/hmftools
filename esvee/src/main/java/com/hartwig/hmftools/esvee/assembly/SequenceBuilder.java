package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.common.utils.Doubles.medianInteger;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REPEAT_MAX_BASE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_LONG_REPEAT_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_MEDIUM_REPEAT_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_MED_QUAL_PENALTY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_LOW_REPEAT_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_LENGTH_1;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_LENGTH_2;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_LENGTH_3;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_PENALTY_1;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_PENALTY_2;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_PENALTY_3;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.READ_MISMATCH_PENALTY_PENALTY_LONG;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.DNA_BASE_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.NO_BASE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.BaseQualType.LOW;
import static com.hartwig.hmftools.esvee.assembly.BaseQualType.MEDIUM;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffInfo.UNSET;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.DELETE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.INSERT;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.MATCH;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHighBaseQual;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public class SequenceBuilder
{
    private final boolean mBuildForwards;
    private final List<ReadParseState> mReads;

    private final boolean mDisableMismatchPenalty;

    private int mCurrentIndex;
    private byte[] mBases;
    private byte[] mBaseQuals;

    private final List<RepeatInfo> mRepeats;

    protected static final int NEXT_BASE_CHECK_COUNT = 3; // needs to be a function of how many skipped indel bases are checked

    public SequenceBuilder(final List<ReadParseState> reads, boolean buildForwards, int initialBaseLength)
    {
        this(reads, buildForwards, initialBaseLength, false);
    }

    public SequenceBuilder(final List<ReadParseState> reads, boolean buildForwards, int initialBaseLength, boolean disableMismatchPenalty)
    {
        mBuildForwards = buildForwards;
        mReads = reads;

        mDisableMismatchPenalty = disableMismatchPenalty;

        mBases = new byte[initialBaseLength];
        mBaseQuals = new byte[initialBaseLength];
        mCurrentIndex = -1;
        mRepeats = Lists.newArrayList();

        buildSequence();
        trimFinalSequence();
        findNonMismatchRepeats();
    }

    public byte[] bases() { return mBases; }
    public int baseLength() { return mBases.length; }
    public String baseString() { return new String(mBases); }
    public byte[] baseQuals() { return mBaseQuals; }
    public List<ReadParseState> reads() { return mReads; }
    public List<RepeatInfo> repeats() { return mRepeats; }

    public RepeatInfo findRepeat(int baseIndex)
    {
        // find a repeat which starts earlier and reaches as far this base
        for(RepeatInfo repeat : mRepeats)
        {
            if(mBuildForwards)
            {
                if(baseIndex >= repeat.Index && baseIndex <= repeat.lastIndex())
                    return repeat;
            }
            else
            {
                if(baseIndex >= repeat.Index && baseIndex <= repeat.lastIndex())
                    return repeat;
            }
        }

        return null;
    }

    private void buildSequence()
    {
        mCurrentIndex = mBuildForwards ? 0 : mBases.length - 1;

        List<ReadParseState> activeReads = Lists.newArrayList(mReads);

        int hpCount = 0;
        byte hpBase = NO_BASE;

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

                byte base = read.currentBase();
                byte qual = read.currentQual();
                boolean isHighQual = isHighBaseQual(qual);

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
                    /*
                    else if(base != consensusBase && belowMinQual(qual))
                    {
                        // low-qual disagreement - ignore regardless of consensus qual
                    }
                    */
                    else
                    {
                        hasMismatch = true;

                        // high-qual mismatch so start tracking frequencies for each base
                        maxQuals = new byte[DNA_BASE_COUNT];
                        highQuals = new int[DNA_BASE_COUNT];
                        mediumQuals = new int[DNA_BASE_COUNT];

                        // back port existing counts to the per-base arrays
                        int baseIndex = AssemblyUtils.baseIndex(consensusBase);
                        maxQuals[baseIndex] = consensusMaxQual;
                        highQuals[baseIndex] = consensusHighQualCount;
                        mediumQuals[baseIndex] = consensusMedQualCount;
                    }
                }

                if(hasMismatch)
                {
                    int baseIndex = AssemblyUtils.baseIndex(base);

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

                boolean tieBreak = false;

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
                    else if(highQuals[b] == maxHighQual && mediumQuals[b] == maxMediumQual)
                    {
                        tieBreak = true;
                    }
                }

                if(tieBreak)
                {
                    // use the read with the lowest mismatches so far
                    double minMismatches = activeReads.stream().mapToDouble(x -> x.mismatchPenalty()).min().orElse(0);

                    ReadParseState read = activeReads.stream().filter(x -> x.mismatchPenalty() == minMismatches).findFirst().orElse(null);

                    if(read != null)
                    {
                        maxBaseIndex = AssemblyUtils.baseIndex(read.currentBase());
                    }
                }

                consensusBase = DNA_BASE_BYTES[maxBaseIndex];
                consensusMaxQual = maxQuals[maxBaseIndex];

                SequenceDiffInfo[] seqDiffInfos = assessDifferences(consensusBase, consensusMaxQual, activeReads, hpBase, hpCount);

                applyReadMismatchesAndConsensus(activeReads, consensusBase, consensusMaxQual, seqDiffInfos);
            }
            else
            {
                mBases[mCurrentIndex] = consensusBase;
                mBaseQuals[mCurrentIndex] = consensusMaxQual;

                mCurrentIndex += mBuildForwards ? 1 : -1;

                for(ReadParseState read : activeReads)
                {
                    read.addBaseMatch(isHighBaseQual(read.currentQual()));
                    read.moveNext();
                }
            }

            if(hpBase == consensusBase)
            {
                ++hpCount;
            }
            else
            {
                hpBase = consensusBase;
                hpCount = 1;
            }

            // purge exhausted reads
            int readIndex = 0;
            while(readIndex < activeReads.size())
            {
                ReadParseState read = activeReads.get(readIndex);

                if(read.exhausted() || read.mismatched())
                    activeReads.remove(readIndex);
                else
                    ++readIndex;
            }

            if(activeReads.isEmpty())
                break;
        }
    }

    private static final int MIN_PREVIOUS_REPEAT_COUNT = READ_MISMATCH_LOW_REPEAT_COUNT - 1;

    private SequenceDiffInfo[] assessDifferences(
            byte consensusBase, byte consensusQual, final List<ReadParseState> activeReads, byte hpBase, int hpCount)
    {
        // assess each read relative to the consensus base to establish if an indel or repeat difference is the cause
        SequenceDiffInfo[] seqDiffInfos = new SequenceDiffInfo[activeReads.size()];

        for(int i = 0; i < seqDiffInfos.length; ++i)
        {
            seqDiffInfos[i] = UNSET;
        }

        int readCount = activeReads.size();

        List<byte[]> allReadNextBases = Lists.newArrayListWithCapacity(readCount);
        List<byte[]> consensusReadNextBases = Lists.newArrayListWithCapacity(readCount);

        SequenceDiffInfo matchInfo = SequenceDiffInfo.fromMatch(-1, mCurrentIndex);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);

            byte[] readNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);

            allReadNextBases.add(readNextBases); // stored even if null

            // no consideration for qual at this point - will revert to a low-qual mismatch after testing indels and repeats
            // boolean baseMatch = basesMatch(consensusBase, read.currentBase(), consensusQual, read.currentQual());

            if(read.currentBase() == consensusBase)
            {
                consensusReadNextBases.add(readNextBases);

                // no need to record read-specific matched info, will not be registered against the read state
                seqDiffInfos[readIndex] = matchInfo;
            }
            else if(readNextBases == null)
            {
                seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(read, mCurrentIndex);
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
                    seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(read, mCurrentIndex);
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
        }

        if(hpCount >= MIN_PREVIOUS_REPEAT_COUNT)
        {
            // look for homopolymer repeat differences
            assessRepeatDifferences(activeReads, seqDiffInfos, hpBase, hpCount);

            if(Arrays.stream(seqDiffInfos).noneMatch(x -> x == UNSET))
                return seqDiffInfos;
        }

        if(consensusNextBases != null)
        {
            assessNovelIndelDifferences(consensusBase, consensusNextBases, activeReads, allReadNextBases, seqDiffInfos);

            if(Arrays.stream(seqDiffInfos).noneMatch(x -> x == UNSET))
                return seqDiffInfos;
        }

        // otherwise check for a repeat expansion / contraction
        assessRepeatDifferences(activeReads, seqDiffInfos, NO_BASE, 0);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            if(seqDiffInfos[readIndex] == UNSET)
            {
                ReadParseState read = activeReads.get(readIndex);
                seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(read, mCurrentIndex);
            }
        }

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

            seqDiffInfos[readIndex] = assessNovelIndelDifferences(
                    read, seqDiffInfo, readNextBases, consensusNextBases, consensusCurrentAndNext);
        }
    }

    protected SequenceDiffInfo assessNovelIndelDifferences(
            final ReadParseState read, final SequenceDiffInfo seqDiffInfo, final byte[] readNextBases,
            final byte[] consensusNextBases, final byte[] consensusCurrentAndNext)
    {
        if(seqDiffInfo.Type == MATCH || seqDiffInfo.Type == BASE || readNextBases == null)
            return seqDiffInfo;

        if(basesMatch(consensusCurrentAndNext, readNextBases))
        {
            // read has inserted the current base
            BaseQualType qualType = read.qualType(read.readIndex());
            String base = String.valueOf((char)read.currentBase());

            return new SequenceDiffInfo(read.readIndex(), mCurrentIndex, base, SequenceDiffType.INSERT, qualType);
        }
        else
        {
            byte[] readCurrentAndNext = currentAndNextBase(read.currentBase(), readNextBases);

            if(basesMatch(consensusNextBases, readCurrentAndNext))
            {
                int previousIndex = read.readIndex() + (mBuildForwards ? -1 : 1);
                BaseQualType qualType = read.straddlingQualType(read.readIndex(), previousIndex); // order doesn't matter

                return new SequenceDiffInfo(
                        read.readIndex(), mCurrentIndex, "", SequenceDiffType.DELETE, qualType);
            }
        }

        return seqDiffInfo;
    }

    private void assessRepeatDifferences(
            final List<ReadParseState> activeReads, final SequenceDiffInfo[] seqDiffInfos, byte hpBase, int hpCount)
    {
        // assume that the repeat ends for at least some of the reads at the base prior to the mismatch base,
        // so search for repeats around this location

        // first of all look for single-base expansions since this is the most common
        String maxRepeatBases = null;
        int maxReadsWithRepeatCount = 0; // number of reads with evidence of the repeat
        int[] maxRepeatReadRepeatCounts = null; // the repeat counts per read for the most likely repeat
        int maxRepeatPreviousRepeatCount = 0;

        if(hpCount > 0)
        {
            maxRepeatBases = String.valueOf((char)hpBase);
            maxRepeatPreviousRepeatCount = hpCount;
            maxRepeatReadRepeatCounts = new int[activeReads.size()];

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);

                int repeatCount = getRepeatCount(read, maxRepeatBases, hpCount);
                maxRepeatReadRepeatCounts[readIndex] = repeatCount;

                if(repeatCount > 0)
                    ++maxReadsWithRepeatCount;
            }
        }
        else
        {
            // search for repeats of 2-base
            for(int repeatLength = 2; repeatLength <= REPEAT_MAX_BASE_LENGTH; ++repeatLength)
            {
                String repeatBases = null;

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

                int repeatIndexStart = mBuildForwards ? mCurrentIndex - 1 : mCurrentIndex + 1;
                int previousRepeatCount = RepeatInfo.getRepeatCount(mBases, repeatBases, repeatIndexStart, !mBuildForwards);

                if(previousRepeatCount < MIN_PREVIOUS_REPEAT_COUNT)
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
        }

        // repeat count must meet minimums to be registered
        int maxObservedRepeatCount = Arrays.stream(maxRepeatReadRepeatCounts).max().orElse(0);

        if(maxObservedRepeatCount < READ_MISMATCH_LOW_REPEAT_COUNT)
            return;

        // must impact at least half the reads, even if exactly 1 of 2 reads - better to assume a repeat-type mismatch
        if(maxReadsWithRepeatCount < 0.5 * activeReads.size())
            return;

        int consensusRepeatCount = findConsensusRepeatCount(maxRepeatReadRepeatCounts);

        // register this consensus repeat
        int repeatBaseLength = maxRepeatBases.length();
        int consensusRepeatOffset = maxRepeatPreviousRepeatCount * repeatBaseLength; // previous repeat bases before the current base
        int consensusRepeatIndexStart = mBuildForwards ? mCurrentIndex - consensusRepeatOffset : mCurrentIndex + consensusRepeatOffset;

        RepeatInfo consensusRepeat = new RepeatInfo(
                consensusRepeatIndexStart, maxRepeatBases, consensusRepeatCount, mBuildForwards ? Orientation.FORWARD : REVERSE);
        mRepeats.add(consensusRepeat);

        // and then register the diff info for each read, including how many repeats it has
        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);
            int readRepeatCount = maxRepeatReadRepeatCounts[readIndex];

            SequenceDiffType diffType = readRepeatCount == consensusRepeatCount ? MATCH : SequenceDiffType.REPEAT;

            int readRepeatStart = SequenceDiffInfo.repeatIndex(
                    read.readIndex(), readRepeatCount, consensusRepeat, consensusRepeatOffset, mBuildForwards, true);

            int readRepeatEnd = SequenceDiffInfo.repeatIndex(
                    read.readIndex(), readRepeatCount, consensusRepeat, consensusRepeatOffset, mBuildForwards, false);

            int readRepeatIndexBegin = mBuildForwards ? readRepeatStart : readRepeatEnd;

            BaseQualType baseQualType = read.rangeMinQualType(readRepeatStart, readRepeatEnd);

            seqDiffInfos[readIndex] = new SequenceDiffInfo(
                    read.readIndex(), mCurrentIndex, maxRepeatBases, diffType, baseQualType, readRepeatCount, readRepeatIndexBegin);
        }
    }

    private int getRepeatCount(final ReadParseState read, final String repeat, int previousRepeatCount)
    {
        return getRepeatCount(read, repeat, previousRepeatCount, mBuildForwards);
    }

    protected static int getRepeatCount(final ReadParseState read, final String repeat, int previousRepeatCount, boolean buildForwards)
    {
        // search forward from the current index
        int countForward = max(RepeatInfo.getRepeatCount(read.read().getBases(), repeat, read.readIndex(), buildForwards), 0);
        return previousRepeatCount + countForward;
    }

    private static int findConsensusRepeatCount(final int[] readRepeatCounts)
    {
        Map<Integer,Integer> repeatFrequencies = Maps.newHashMap();
        List<Integer> allCounts = Lists.newArrayList();

        for(int readIndex = 0; readIndex < readRepeatCounts.length; ++readIndex)
        {
            int readRepeatCount = readRepeatCounts[readIndex];

            if(readRepeatCount < READ_MISMATCH_LOW_REPEAT_COUNT)
                continue;

            Integer readCount = repeatFrequencies.get(readRepeatCount);
            repeatFrequencies.put(readRepeatCount, readCount == null ? 1 : readCount + 1);
            allCounts.add(readRepeatCount);
        }

        int maxRepeatFreq = 0;
        int maxRepeatCount = 0;

        // favour a shorter repeat count if a tie break? or if just 2 reads?
        for(Map.Entry<Integer,Integer> entry : repeatFrequencies.entrySet())
        {
            if(entry.getValue() > maxRepeatFreq || (entry.getValue() == maxRepeatFreq && entry.getKey() < maxRepeatCount))
            {
                maxRepeatCount = entry.getKey();
                maxRepeatFreq = entry.getValue();
            }
        }

        if(maxRepeatFreq < 0.33 * allCounts.size())
        {
            return medianInteger(allCounts);
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

        // only take the next consensus bases if it is present in the majority of consensus reads
        double matchedFraction = consensusReadNextBaseMatchCount / (double)consensusReadCount;

        return matchedFraction > 0.5 ? consensusNextBases : null;
    }

    private void applyReadMismatchesAndConsensus(
            final List<ReadParseState> activeReads, byte consensusBase, byte consensusQual, final SequenceDiffInfo[] seqDiffInfos)
    {
        // the last repeat added corresponds to the any read(s) with a repeat mismatch
        RepeatInfo consensusRepeat = Arrays.stream(seqDiffInfos).anyMatch(x -> x.RepeatCount > 0) ? mRepeats.get(mRepeats.size() - 1) : null;

        int previousReadBases = consensusRepeat != null ? abs(mCurrentIndex - consensusRepeat.Index) : 0;

        setConsensusFromMismatches(consensusBase, consensusQual, consensusRepeat);

        // move each read's current index onwards or pause as required by differences to the consensus
        // also register mismatch penalties based on type, seq-tech, quals, repeat-diffs etc
        int readIndex = 0;
        while(readIndex < activeReads.size())
        {
            ReadParseState read = activeReads.get(readIndex);
            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];

            applyReadMismatches(read, seqDiffInfo, consensusRepeat, previousReadBases);

            if(read.mismatched())
            {
                activeReads.remove(readIndex);
            }
            else
            {
                ++readIndex;
            }
        }
    }

    private void setConsensusFromMismatches(byte consensusBase, byte consensusQual, final RepeatInfo consensusRepeat)
    {
        if(consensusRepeat == null)
        {
            mBases[mCurrentIndex] = consensusBase;
            mBaseQuals[mCurrentIndex] = consensusQual;
            mCurrentIndex += mBuildForwards ? 1 : -1;
            return;
        }

        int repeatLength = consensusRepeat.repeatLength();

        // adjust each read's index relative to the consensus repeat
        int previousReadBases = abs(mCurrentIndex - consensusRepeat.Index);
        int remainingRepeatBases = consensusRepeat.totalLength() - previousReadBases;

        int repeatBaseIndex = mBuildForwards ? 0 : repeatLength - 1;
        String repeatBases = consensusRepeat.Bases;

        // apply the repeat to the consensus bases
        for(int i = 0; i < remainingRepeatBases; ++i)
        {
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
    }

    protected void applyReadMismatches(
            final ReadParseState read, final SequenceDiffInfo seqDiffInfo, final RepeatInfo consensusRepeat, final int previousRepeatLength)
    {
        // register mismatches and move the read index on as required by the match type
        if(seqDiffInfo.Type != MATCH)
        {
            seqDiffInfo.MismatchPenalty = calcMismatchPenalty(seqDiffInfo, consensusRepeat);
            read.addMismatchInfo(seqDiffInfo);

            if(checkReadMismatches(read))
            {
                read.markMismatched();
                return;
            }
        }
        else
        {
            read.addBaseMatch(isHighBaseQual(read.currentQual()));
        }

        read.moveOnMatchType(seqDiffInfo);

        if(consensusRepeat != null)
        {
            processReadRepeat(read, seqDiffInfo, consensusRepeat, previousRepeatLength);
        }
    }

    protected void processReadRepeat(
            final ReadParseState read, final SequenceDiffInfo seqDiffInfo, final RepeatInfo consensusRepeat, final int previousRepeatLength)
    {
        int repeatedMatches = consensusRepeat.totalLength();

        // move each read's index based on diffs vs the consensus repeat
        int remainingReadRepeatBases = seqDiffInfo.RepeatCount * consensusRepeat.repeatLength() - previousRepeatLength;
        read.moveNextBases(remainingReadRepeatBases);

        int readRepeatStart = seqDiffInfo.repeatIndex(consensusRepeat, mBuildForwards, true);
        int readRepeatEnd = seqDiffInfo.repeatIndex(consensusRepeat, mBuildForwards, false);

        int highQualRepeatBases = read.highQualBaseCount(readRepeatStart, readRepeatEnd);
        read.addBaseMatches(repeatedMatches, highQualRepeatBases);
    }

    private boolean checkReadMismatches(final ReadParseState read)
    {
        if(mDisableMismatchPenalty)
            return false;

        // CHECK and compare with AssemblyUtils.mismatchesPerComparisonLength
        return exceedsReadMismatches(read);
    }

    public static boolean exceedsReadMismatches(final ReadParseState read)
    {
        double mismatches = read.mismatchPenalty();
        double permittedPenalty = permittedReadMismatches(read.overlapBaseCount());
        return mismatches > permittedPenalty;
    }

    public static double permittedReadMismatches(int readOverlap)
    {
        if(readOverlap <= READ_MISMATCH_PENALTY_LENGTH_1)
            return READ_MISMATCH_PENALTY_PENALTY_1;

        if(readOverlap <= READ_MISMATCH_PENALTY_LENGTH_2)
            return READ_MISMATCH_PENALTY_PENALTY_2;

        if(readOverlap <= READ_MISMATCH_PENALTY_LENGTH_3)
            return READ_MISMATCH_PENALTY_PENALTY_3;

        return READ_MISMATCH_PENALTY_PENALTY_LONG;
    }

    public static double calcMismatchPenalty(final SequenceDiffInfo seqDiffInfo, @Nullable final RepeatInfo consensusRepeat)
    {
        if(seqDiffInfo.Type == MATCH)
            return 0;

        if(seqDiffInfo.QualType == LOW)
            return 0;

        if(seqDiffInfo.Type == BASE || seqDiffInfo.Type == INSERT || seqDiffInfo.Type == DELETE)
        {

            if(seqDiffInfo.QualType == MEDIUM)
                return READ_MISMATCH_MED_QUAL_PENALTY;
            else
                return READ_MISMATCH_PENALTY;
        }

        // repeat adjustments
        if(consensusRepeat == null)
            return READ_MISMATCH_PENALTY;

        int repeatCount = consensusRepeat.Count;
        int repeatCountDiff = abs(seqDiffInfo.RepeatCount - consensusRepeat.Count);

        int penaltyCount;
        if(repeatCount < READ_MISMATCH_MEDIUM_REPEAT_COUNT)
        {
            penaltyCount = max(repeatCountDiff - 1, 0);
        }
        else if(repeatCount < READ_MISMATCH_LONG_REPEAT_COUNT)
        {
            penaltyCount = max(repeatCountDiff - 2, 0);
        }
        else
        {
            penaltyCount = 0;
        }

        if(seqDiffInfo.QualType == BaseQualType.MEDIUM)
            return penaltyCount * READ_MISMATCH_MED_QUAL_PENALTY;

        return penaltyCount;
    }

    public void populateBases(final byte[] bases, final int indexStart)
    {
        int baseIndex = indexStart;
        for(int i = 0; i < bases.length; ++i)
        {
            int destIndex = mBuildForwards ? i : bases.length - i - 1;
            bases[destIndex] = mBases[baseIndex];

            baseIndex += mBuildForwards ? 1 : -1;

            if(baseIndex < 0 || baseIndex >= mBases.length)
                break;
        }
    }

    private void findNonMismatchRepeats()
    {
        // find repeats other than those resulting from read mismatches
        List<RepeatInfo> allRepeats = RepeatInfo.findRepeats(mBases);

        if(allRepeats.isEmpty())
            return;

        // keep any which do not overlap with the mismatch repeats, and favour longer counts
        Collections.sort(allRepeats, Comparator.comparingInt(x -> -x.Count));
        int index = 0;

        while(index < allRepeats.size())
        {
            RepeatInfo repeat = allRepeats.get(index);

            if(mRepeats.stream().anyMatch(x -> x.overlaps(repeat)))
            {
                allRepeats.remove(index);
            }
            else
            {
                if(mBuildForwards)
                {
                    mRepeats.add(repeat);
                }
                else
                {
                    mRepeats.add(new RepeatInfo(repeat.lastIndex(), repeat.Bases, repeat.Count, REVERSE));
                }
            }
        }
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

        if(!mBuildForwards)
        {
            // adjust the repead indicies down
            for(int i = 0; i < mRepeats.size(); ++i)
            {
                RepeatInfo repeat = mRepeats.get(i);
                mRepeats.set(i, new RepeatInfo(repeat.Index - reduction, repeat.Bases, repeat.Count, repeat.direction()));
            }
        }
    }
}
