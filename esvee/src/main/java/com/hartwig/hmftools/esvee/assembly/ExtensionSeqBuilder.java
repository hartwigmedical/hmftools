package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.N_BASE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.compareSequences;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.getReadIndexAtReferencePosition;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class ExtensionSeqBuilder
{
    private final Junction mJunction;
    private final List<ReadState> mReads;
    private final boolean mIsForward;
    private final int mMaxMismatches;

    private byte[] mBases;
    private byte[] mBaseQuals;
    private int[] mReadCount;
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

            if(readJunctionIndex == INVALID_INDEX)
                continue;

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
        mMinSupportLength = 0;
        mExtensionRepeats = Lists.newArrayList();

        mIsValid = true;

        buildConsensus();

        assignReads();

        determineFinalBases();
    }

    public byte[] extensionBases() { return mBases; }
    public byte[] baseQualitiies() { return mBaseQuals; }
    public int minSupportLength() { return mMinSupportLength; }
    public List<RepeatInfo> repeatInfo() { return mExtensionRepeats; }
    public boolean isValid() { return mIsValid; }

    public List<SupportRead> formAssemblySupport()
    {
        List<SupportRead> supportReads = Lists.newArrayList();

        for(ReadState read : mReads)
        {
            int permittedReadMismatches = mismatchesPerComparisonLength(read.extensionLength());

            if(read.Mismatches > permittedReadMismatches)
                continue;

            supportReads.add(new SupportRead(
                    read.read(), SupportType.JUNCTION, read.junctionIndex(), read.matchedBases(), read.Mismatches));
        }

        return supportReads;
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

                if(base == N_BASE)
                {
                    base = DNA_BASE_BYTES[0];
                    qual = 0;
                }

                read.moveNext(mIsForward);

                if(readCounts == null)
                {
                    if(consensusBase == 0
                    || AssemblyUtils.basesMatch(base, consensusBase, (byte)qual, (byte)consensusMaxQual, LOW_BASE_QUAL_THRESHOLD))
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
                consensusReadCount = readCounts[maxBaseIndex];
            }

            mBases[extensionIndex] = consensusBase;
            mBaseQuals[extensionIndex] = (byte)consensusMaxQual;
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
                byte qual = read.currentQual();

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

    private void determineFinalBases()
    {
        int maxValidExtensionLength = 0;

        for(ReadState read : mReads)
        {
            if(read.Mismatches > mMaxMismatches)
                continue;

            // note the read's extension length does not include the junction base itself, hence the +1
            maxValidExtensionLength = max(read.extensionLength() + 1, maxValidExtensionLength);
        }

        if(maxValidExtensionLength == 0)
        {
            mIsValid = false;
            return;
        }

        // trim extension bases if required
        if(maxValidExtensionLength == mBases.length)
            return;

        int reduction = mBases.length - maxValidExtensionLength;
        int startIndex = mIsForward ? 0 : reduction;
        int endIndex = mIsForward ? maxValidExtensionLength - 1 : mBases.length - 1;
        mBases = subsetArray(mBases, startIndex, endIndex);
        mBaseQuals = subsetArray(mBaseQuals, startIndex, endIndex);

        mExtensionRepeats.clear();
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);

        if(repeats != null)
            mExtensionRepeats.addAll(repeats);
    }

    public static int calcReadSequenceMismatches(
            final boolean isForward, final byte[] extensionBases, final byte[] extensionQuals, final List<RepeatInfo> extensionRepeats,
            final Read read, final int readJunctionIndex, final int maxMismatches)
    {
        int readStartIndex = isForward ? readJunctionIndex : 0;
        int readEndIndex = isForward ? read.basesLength() - 1 : readJunctionIndex;

        // for -ve orientations, if extension sequence length = 10, with 0-8 being soft-clip and 9 being the first ref and junction index
        // and the read has 5 bases of soft-clip then read's start index will be 0 -> 4 + 1 = 5
        // so the comparison offset in the extension sequence is
        int extSeqReadStartIndex = isForward ? 0 : extensionBases.length - 1 - readJunctionIndex;

        byte[] readExtensionBases = subsetArray(read.getBases(), readStartIndex, readEndIndex);
        byte[] readExtensionQuals = subsetArray(read.getBaseQuality(), readStartIndex, readEndIndex);
        List<RepeatInfo> readRepeats = RepeatInfo.findRepeats(readExtensionBases);

        return compareSequences(
                extensionBases, extensionQuals, extSeqReadStartIndex, extensionBases.length - 1, extensionRepeats,
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
        public int extensionLength() { return mExtensionLength; }
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
                    mRead.id(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                    mExtensionLength, mCurrentIndex, Mismatches, mExhausted ? "exhausted" : "");
        }
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mBases); }
}