package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_REPEAT_BASE_COUNT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.DNA_BASE_COUNT;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHighBaseQual;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;

public class SequenceBuilder
{
    private final boolean mBuildForwards;
    private final int mInitialBaseLength;
    private final List<ReadParseState> mReads;

    private byte[] mBases;
    private byte[] mBaseQuals;

    public SequenceBuilder(final List<ReadParseState> reads, boolean buildForwards, int initialBaseLength)
    {
        mBuildForwards = buildForwards;
        mInitialBaseLength = initialBaseLength;
        mReads = reads;

        mBases = new byte[initialBaseLength];
        mBaseQuals = new byte[initialBaseLength];

        buildSequence();
    }

    private void buildSequence()
    {
        int extensionIndex = mBuildForwards ? 0 : mBases.length - 1;

        List<ReadParseState> activeReads = Lists.newArrayList(mReads);

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            byte consensusBase = 0;
            byte consensusMaxQual = 0;
            int consensusHighQualCount = 0;
            int consensusMedQualCount = 0;

            byte[] maxQuals = null; // max qual per base
            int[] highQuals = null; // count of reads with high-quals per base
            int[] mediumQuals = null; // count of reads with high-quals per base

            boolean hasActiveReads = false;

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

                SequenceDiffInfo[] seqDiffInfos = assessDifferences(extensionIndex, consensusBase, activeReads);
            }

            mBases[extensionIndex] = consensusBase;
            mBaseQuals[extensionIndex] = consensusMaxQual;

            // progress the read index and purge exhausted reads
            int readIndex = 0;
            while(readIndex < activeReads.size())
            {
                ReadParseState read = activeReads.get(readIndex);
                read.moveNext();

                if(read.exhausted())
                    activeReads.remove(readIndex);
                else
                    ++readIndex;
            }

            if(activeReads.isEmpty())
                break;

            if(mBuildForwards)
                ++extensionIndex;
            else
                --extensionIndex;
        }
    }

    private static final int NEXT_BASE_CHECK_COUNT = 3;

    private SequenceDiffInfo[] assessDifferences(int extensionIndex, byte consensusBase, final List<ReadParseState> activeReads)
    {
        // assess each read relative to the consensus base to establish if an indel or repeat difference is the cause
        SequenceDiffInfo[] seqDiffInfos = new SequenceDiffInfo[activeReads.size()];

        int readCount = activeReads.size();

        // first check if bases match after this single base mismatch, ie no indel-based differences
        List<byte[]> allReadNextBases = Lists.newArrayListWithCapacity(readCount);
        List<byte[]> consensusReadNextBases = Lists.newArrayListWithCapacity(readCount);

        for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
        {
            ReadParseState read = activeReads.get(readIndex);

            byte[] readNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);

            allReadNextBases.add(readNextBases); // stored even if null

            if(read.currentBase() == consensusBase)
                consensusReadNextBases.add(readNextBases);
        }

        // first check if the consensus base matches the prior base for the reads which match
        byte[] consensusNextBases = findConsensusNextBases(consensusReadNextBases);

        SequenceDiffInfo matchInfo = SequenceDiffInfo.fromMatch(extensionIndex);

        int consensusMatchCount = 0;

        if(consensusNextBases != null)
        {
            int postConsensusMismatchCount = 0;

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);
                byte[] readNextBases = allReadNextBases.get(readIndex);

                if(readNextBases == null)
                {
                    if(consensusBase == read.currentBase())
                    {
                        ++consensusMatchCount;
                        seqDiffInfos[readIndex] = matchInfo;
                    }

                    continue;
                }

                if(!basesMatch(readNextBases, consensusNextBases))
                {
                    ++postConsensusMismatchCount;
                }
                else if(consensusBase != read.currentBase() && aboveMinQual(read.currentQual()))
                {
                    seqDiffInfos[readIndex] = SequenceDiffInfo.fromSnv(extensionIndex, read.currentBase());
                }
                else
                {
                    ++consensusMatchCount;
                    seqDiffInfos[readIndex] = matchInfo;
                }
            }

            if(postConsensusMismatchCount == 0)
            {
                // difference is explained by an SNV in some of the reads
                return seqDiffInfos;
            }

            assessNovelIndelDifferences(extensionIndex, consensusBase, consensusNextBases, activeReads, allReadNextBases, seqDiffInfos);

            int remainingMismatchCount = (int) Arrays.stream(seqDiffInfos).filter(x -> x == null).count();

            if(remainingMismatchCount == 0)
                return seqDiffInfos;
        }

        // otherwise check for a repeat expansion / contraction
        assessRepeatDifferences(extensionIndex, activeReads, seqDiffInfos);

        return seqDiffInfos;
    }

    private void assessNovelIndelDifferences(
            int extensionIndex, byte consensusBase, final byte[] consensusNextBases, final List<ReadParseState> activeReads,
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

            if(seqDiffInfo.Type == SequenceDiffType.MATCH || readNextBases == null)
                continue;

            if(basesMatch(consensusCurrentAndNext, readNextBases))
            {
                seqDiffInfos[readIndex] = new SequenceDiffInfo(extensionIndex, null, SequenceDiffType.DELETE, 0);
            }
            else
            {
                byte[] readCurrentAndNext = currentAndNextBase(read.currentBase(), readNextBases);
                if(basesMatch(consensusNextBases, readCurrentAndNext))
                {
                    seqDiffInfos[readIndex] = new SequenceDiffInfo(
                            extensionIndex, String.valueOf((char)read.currentBase()), SequenceDiffType.INSERT, 0);
                }
            }
        }
    }

    private void assessRepeatDifferences(int extensionIndex, final List<ReadParseState> activeReads, final SequenceDiffInfo[] seqDiffInfos)
    {
        // assume that the repeat ends for at least some of the reads at the base prior to the mismatch base,
        // so search for repeats around this location

        // first of all look for single-base expansions since this is the most common
        String maxRepeatBases = null;
        int maxReadRepeatCount = 0; // number of reads with evidence of the repeat
        int[] maxRepeatReadRepeatCounts = null; // the repeat counts per read for the most likely repeat

        for(int repeatLength = 1; repeatLength <= MAX_REPEAT_BASE_COUNT; ++repeatLength)
        {
            String repeatBases = null;
            int[] repeatReadRepeatCounts = new int[activeReads.size()];
            int nonZeroRepeatReadCount = 0;

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);

                if(readIndex == 0)
                {
                    if(repeatLength == 1)
                    {
                        byte previousBase = read.getPreviousBase();
                        repeatBases = String.valueOf((char)previousBase);
                    }
                    else
                    {
                        byte[] previousBases = read.getPreviousBases(repeatLength);
                        repeatBases = new String(previousBases);
                    }
                }

                // find repeat count for this read
                int repeatIndexStart = mBuildForwards ? read.readIndex() - repeatBases.length() : read.readIndex() + repeatBases.length();

                int repeatCount = getRepeatCount(read.read(), repeatBases, repeatIndexStart);
                repeatReadRepeatCounts[readIndex] = repeatCount;

                if(repeatCount > 0)
                    ++nonZeroRepeatReadCount;
            }

            if(nonZeroRepeatReadCount > maxReadRepeatCount)
            {
                maxReadRepeatCount = nonZeroRepeatReadCount;
                maxRepeatReadRepeatCounts = repeatReadRepeatCounts;
                maxRepeatBases = repeatBases;
            }
        }

        double readsWithRepeatFraction = maxReadRepeatCount / (double)activeReads.size();

        if(readsWithRepeatFraction > 0.5)
        {
            int consensusRepeatCount = findConsensusRepeatCount(seqDiffInfos, maxRepeatReadRepeatCounts);

            for(int readIndex = 0; readIndex < activeReads.size(); ++readIndex)
            {
                ReadParseState read = activeReads.get(readIndex);
                SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];
                int readRepeatCount = maxRepeatReadRepeatCounts[readIndex];

                SequenceDiffType diffType = seqDiffInfo != null && seqDiffInfo.Type == SequenceDiffType.MATCH ?
                        SequenceDiffType.MATCH : SequenceDiffType.REPEAT;

                int readRepeatDiff = readRepeatCount - consensusRepeatCount;

                seqDiffInfos[readIndex] = new SequenceDiffInfo(extensionIndex, maxRepeatBases, diffType, readRepeatDiff);
            }
        }
    }

    private static int getRepeatCount(final Read read, final String repeat, int readIndexAnchor)
    {
        int countBack = RepeatInfo.getRepeatCount(read, repeat, readIndexAnchor, false);
        int countForward = RepeatInfo.getRepeatCount(read, repeat, readIndexAnchor + repeat.length(), true);
        return countBack + countForward;
    }

    private static int findConsensusRepeatCount(final SequenceDiffInfo[] seqDiffInfos, final int[] readRepeatCounts)
    {
        Map<Integer,Integer> repeatFrequencies = Maps.newHashMap();

        for(int readIndex = 0; readIndex < readRepeatCounts.length; ++readIndex)
        {
            SequenceDiffInfo seqDiffInfo = seqDiffInfos[readIndex];

            if(seqDiffInfo.Type == SequenceDiffType.MATCH)
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

    private static boolean basesMatch(final byte[] bases1, final byte[] bases2)
    {
        if(bases1.length != bases2.length)
            return false;

        for(int i = 0; i < bases1.length; ++i)
        {
            if(bases1[i] != bases2[i])
                return false;
        }

        return true;
    }
}
