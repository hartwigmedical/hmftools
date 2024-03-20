package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.compareSequences;
import static com.hartwig.hmftools.esvee.read.ReadUtils.getReadIndexAtReferencePosition;
import static com.hartwig.hmftools.esvee.read.ReadUtils.subsetArray;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.RepeatInfo;
import com.hartwig.hmftools.esvee.common.SupportType;
import com.hartwig.hmftools.esvee.read.Read;

public class ExtensionSeqBuilder
{
    private final Junction mJunction;
    private final List<ReadState> mReads;
    private final boolean mIsForward;
    private final int mMaxMismatches;

    private byte[] mBases;
    private byte[] mBaseQuals;
    private int[] mReadCount;
    private int[] mBaseQualTotals;
    private int mMinSupportLength;

    private List<RepeatInfo> mExtensionRepeats;

    private boolean mIsValid;

    public ExtensionSeqBuilder(final Junction junction, final List<Read> reads, final int maxMismatches)
    {
        mJunction = junction;
        mIsForward = mJunction.isForward();
        mMaxMismatches = maxMismatches;
        mReads = Lists.newArrayListWithCapacity(reads.size());

        int maxExtension = 0;

        for(Read read : reads)
        {
            int readJunctionIndex = getReadIndexAtReferencePosition(read, junction.Position, true);

            // calculate how many bases beyond the junction the read extends
            // for positive orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 3
            // for negative orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            int extensionLength = junction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

            maxExtension = max(maxExtension, extensionLength);

            mReads.add(new ReadState(read, readJunctionIndex, extensionLength));
        }

        int baseLength = maxExtension + 1;

        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];
        mReadCount = new int[baseLength];
        mBaseQualTotals = new int[baseLength];
        mMinSupportLength = 0;
        mExtensionRepeats = Lists.newArrayList();

        mIsValid = true;

        buildConsensus();

        assignReads();
    }

    public byte[] extensionBases() { return mBases; }
    public byte[] baseQualitiies() { return mBaseQuals; }
    public int[] baseQualTotals() { return mBaseQualTotals; }
    public int minSupportLength() { return mMinSupportLength; }
    public List<RepeatInfo> repeatInfo() { return mExtensionRepeats; }
    public boolean isValid() { return mIsValid; }

    public List<AssemblySupport> formAssemblySupport()
    {
        return mReads.stream().filter(x -> x.Mismatches <= mMaxMismatches)
                .map(x -> new AssemblySupport(
                        x.read(), SupportType.JUNCTION, 0, x.junctionIndex(), x.matchedBases(), x.Mismatches))
                .collect(Collectors.toList());
    }

    public int mismatches() { return (int)mReads.stream().filter(x -> x.Mismatches > mMaxMismatches).count(); }

    private void buildConsensus()
    {
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        int baseCount = Nucleotides.DNA_BASES.length;

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            byte consensusBase = 0;
            int consensusMaxQual = 0;
            int consensusQualTotal = 0;
            int consensusReadCount = 0;

            // per-base arrays are only used for high-qual mismatches
            int[] readCounts = null;
            int[] totalQuals = null;
            int[] maxQuals = null;

            for(ReadState read : mReads)
            {
                if(read.exhausted())
                    continue;

                byte base = read.currentBase();
                int qual = read.currentQual();

                read.moveNext(mIsForward);

                if(readCounts == null)
                {
                    if(consensusBase == 0 || base == consensusBase
                    || qual < LOW_BASE_QUAL_THRESHOLD || consensusMaxQual < LOW_BASE_QUAL_THRESHOLD)
                    {
                        if(consensusBase == 0)
                        {
                            consensusBase = base;
                            consensusMaxQual = qual;
                        }
                        else
                        {
                            consensusMaxQual = max(qual, consensusMaxQual);
                        }

                        consensusQualTotal += qual;
                        ++consensusReadCount;
                        continue;
                    }

                    // start tracking frequencies for each base
                    readCounts = new int[baseCount];
                    totalQuals = new int[baseCount];
                    maxQuals = new int[baseCount];

                    // back port existing counts to the per-base arrays
                    int baseIndex = Nucleotides.baseIndex(consensusBase);
                    totalQuals[baseIndex] = consensusQualTotal;
                    maxQuals[baseIndex] = consensusMaxQual;
                    readCounts[baseIndex] = consensusReadCount;
                }

                int baseIndex = Nucleotides.baseIndex(base);

                totalQuals[baseIndex] += qual;
                maxQuals[baseIndex] = max(maxQuals[baseIndex], qual);
                ++readCounts[baseIndex];
            }

            if(readCounts != null)
            {
                // take the bases with the highest qual totals, and stop if the read count falls below the min
                int maxQual = 0;
                int maxBaseIndex = 0;
                for(int b = 0; b < baseCount; ++b)
                {
                    if(totalQuals[b] > maxQual)
                    {
                        maxQual = totalQuals[b];
                        maxBaseIndex = b;
                    }
                }

                consensusBase = DNA_BASE_BYTES[maxBaseIndex];
                consensusMaxQual = (byte)maxQuals[maxBaseIndex];
                consensusQualTotal = totalQuals[maxBaseIndex];
                consensusReadCount = readCounts[maxBaseIndex];
            }

            mBases[extensionIndex] = consensusBase;
            mBaseQuals[extensionIndex] = (byte)consensusMaxQual;
            mBaseQualTotals[extensionIndex] = consensusQualTotal;
            mReadCount[extensionIndex] = consensusReadCount;

            if(consensusReadCount >= PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                ++mMinSupportLength;

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;
        }

        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);

        if(repeats != null)
            mExtensionRepeats.addAll(repeats);
    }

    private void assignReads()
    {
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        mReads.forEach(x -> x.resetIndex());

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            for(ReadState read : mReads)
            {
                if(read.exhausted())
                    continue;

                byte base = read.currentBase();
                int qual = read.currentQual();

                if(qual >= LOW_BASE_QUAL_THRESHOLD && mBaseQuals[extensionIndex] >= LOW_BASE_QUAL_THRESHOLD)
                {
                    if(base != mBases[extensionIndex])
                        ++read.Mismatches;
                }

                read.moveNext(mIsForward);
            }

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;
        }

        List<ReadState> mismatchReads = mReads.stream().filter(x -> x.Mismatches > mMaxMismatches).collect(Collectors.toList());

        if(mismatchReads.isEmpty())
            return;

        // test mismatch reads using the sequence matcher
        mismatchReads.forEach(x -> x.resetIndex());
        mismatchReads.forEach(x -> x.Mismatches = 0);

        for(ReadState read : mismatchReads)
        {
            read.Mismatches = calcReadSequenceMismatches(
                    mIsForward, mBases, mBaseQuals, mExtensionRepeats, read.read(), read.junctionIndex(), mMaxMismatches);
        }
    }

    public static int calcReadSequenceMismatches(
            final boolean isForward, final byte[] extensionBases, final byte[] extensionQuals, final List<RepeatInfo> extensionRepeats,
            final Read read, final int readJunctionIndex, final int maxMismatches)
    {
        int readStartIndex = isForward ? readJunctionIndex : 0;
        int readEndIndex = isForward ? read.basesLength() - 1 : readJunctionIndex;
        byte[] readExtensionBases = subsetArray(read.getBases(), readStartIndex, readEndIndex);
        byte[] readExtensionQuals = subsetArray(read.getBaseQuality(), readStartIndex, readEndIndex);
        List<RepeatInfo> readRepeats = RepeatInfo.findRepeats(readExtensionBases);

        return compareSequences(
                extensionBases, extensionQuals, 0, extensionBases.length - 1, extensionRepeats,
                readExtensionBases, readExtensionQuals, 0, readExtensionBases.length - 1,
                readRepeats != null ? readRepeats : Collections.emptyList(), maxMismatches);
    }

    private class ReadState
    {
        private final Read mRead;
        private final int mJunctionIndex;
        private final int mExtensionLength;
        private int mCurrentIndex;
        private boolean mExhausted;

        public int Mismatches;

        public ReadState(final Read read, final int junctionIndex, final int extensionLength)
        {
            mRead = read;
            mJunctionIndex = junctionIndex;
            mCurrentIndex = junctionIndex;
            mExtensionLength = extensionLength;
            mExhausted = false;
            Mismatches = 0;
        }

        public Read read() { return mRead; }
        public int junctionIndex() { return mJunctionIndex; }
        public int matchedBases() { return mExtensionLength - Mismatches; }
        public boolean exhausted() { return mExhausted; }

        public byte currentBase() { return mRead.getBases()[mCurrentIndex]; }
        public byte currentQual() { return mRead.getBaseQuality()[mCurrentIndex]; }

        public void moveNext(boolean isForward)
        {
            if(isForward)
                ++mCurrentIndex;
            else
                --mCurrentIndex;

            mExhausted = mCurrentIndex < 0 || mCurrentIndex >= mRead.basesLength();
        }

        public void resetIndex()
        {
            mCurrentIndex = mJunctionIndex;
            mExhausted = false;
        }

        public String toString()
        {
            return format("%s: range(%d - %d) cigar(%s) extLen(%d) curIndex(%d) mismatches(%d) %s",
                    mRead.getName(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                    mExtensionLength, mCurrentIndex, Mismatches, mExhausted ? "exhausted" : "");
        }
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mBases); }
}