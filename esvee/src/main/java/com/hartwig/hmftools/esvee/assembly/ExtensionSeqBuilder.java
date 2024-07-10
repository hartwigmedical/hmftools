package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.N_BASE;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.getReadIndexAtReferencePosition;

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

    private byte[] mBases;
    private byte[] mBaseQuals;
    private int mMinSupportLength;

    private List<RepeatInfo> mExtensionRepeats;

    private boolean mIsValid;

    public ExtensionSeqBuilder(final Junction junction, final List<Read> reads)
    {
        mJunction = junction;
        mIsForward = mJunction.isForward();
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

            // indel-based reads may have much shorter extensions than the typical 32-base requirement, so ensure they can extend far enough
            // to satisfy the min-high-qual match filter used post-consensus
            if(extensionLength < PRIMARY_ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                continue;

            maxExtension = max(maxExtension, extensionLength);

            mReads.add(new ReadState(read, readJunctionIndex, extensionLength));
        }

        int baseLength = maxExtension + 1; // since the junction base itself is included (which is a ref base)

        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];
        mMinSupportLength = 0;
        mExtensionRepeats = Lists.newArrayList();

        mIsValid = true;

        buildSequence();

        formConsensusSequence();

        findRepeats();

        assignReads();

        determineFinalBases();
    }

    public byte[] extensionBases() { return mBases; }
    public byte[] baseQualities() { return mBaseQuals; }
    public int minSupportLength() { return mMinSupportLength; }
    public List<RepeatInfo> repeatInfo() { return mExtensionRepeats; }
    public boolean isValid() { return mIsValid; }

    public List<SupportRead> formAssemblySupport()
    {
        List<SupportRead> supportReads = Lists.newArrayList();

        for(ReadState read : mReads)
        {
            if(read.exceedsMaxMismatches())
                continue;

            if(read.highQualMatches() < PRIMARY_ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                continue;

            supportReads.add(new SupportRead(
                    read.read(), SupportType.JUNCTION, read.junctionIndex(), read.matchedBases(), read.mismatches()));
        }

        return supportReads;
    }

    public List<Read> mismatchReads()
    {
        return mReads.stream().filter(x -> x.exceedsMaxMismatches()).map(x -> x.mRead).collect(Collectors.toList());
    }

    public int mismatches() { return (int)mReads.stream().filter(x -> x.exceedsMaxMismatches()).count(); }

    private void buildSequence()
    {
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        int baseCount = Nucleotides.DNA_BASES.length;

        mMinSupportLength = 0;

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

                if(read.exceedsMaxMismatches()) // be relevant when building the final consensus
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

            if(consensusReadCount >= PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
                ++mMinSupportLength;

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;
        }
    }

    private void formConsensusSequence()
    {
        // test reads against the initial sequence, filtering out those with too many mismatches, and forming a consensus from the rest
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
                        read.addMismatch();
                    else
                        read.addHighQualMatch();
                }

                read.moveNext(mIsForward);
            }

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;
        }

        int mismatchReadCount = (int)mReads.stream().filter(x -> x.exceedsMaxMismatches()).count();

        // if no reads have been excluded then keep the initial sequence
        if(mismatchReadCount == 0)
            return;

        // otherwise filter out the mismatch reads and build the sequence again
        int baseLength = mBases.length;
        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];

        mReads.forEach(x -> x.resetIndex());

        buildSequence();
    }

    private void findRepeats()
    {
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);

        if(repeats != null)
            mExtensionRepeats.addAll(repeats);
    }

    private void assignReads()
    {
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        for(ReadState read : mReads)
        {
            read.resetIndex();

            if(!read.exceedsMaxMismatches())
                read.resetMatches();
        }

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            for(ReadState read : mReads)
            {
                if(read.exhausted() || read.exceedsMaxMismatches())
                    continue;

                byte base = read.currentBase();
                byte qual = read.currentQual();

                if(qual >= LOW_BASE_QUAL_THRESHOLD && mBaseQuals[extensionIndex] >= LOW_BASE_QUAL_THRESHOLD)
                {
                    if(base != mBases[extensionIndex])
                        read.addMismatch();
                    else
                        read.addHighQualMatch();
                }

                read.moveNext(mIsForward);
            }

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;
        }

        // no longer use the sequence matcher on reads exceeding the mismatch limit - this is just handled when re-testing all junction reads
    }

    private void determineFinalBases()
    {
        int maxValidExtensionLength = 0;

        for(ReadState read : mReads)
        {
            if(read.exceedsMaxMismatches())
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

    private class ReadState
    {
        private final Read mRead;
        private final int mJunctionIndex;
        private final int mExtensionLength;
        private int mCurrentIndex;
        private boolean mExhausted;
        private final int mPermittedMismatches;

        private int mMismatches;
        public int mHighQualMatches;

        public ReadState(final Read read, final int junctionIndex, final int extensionLength)
        {
            mRead = read;
            mJunctionIndex = junctionIndex;
            mCurrentIndex = junctionIndex;
            mExtensionLength = extensionLength;
            mPermittedMismatches = mismatchesPerComparisonLength(extensionLength);
            mExhausted = false;

            mMismatches = 0;
            mHighQualMatches = 0;
        }

        public Read read() { return mRead; }
        public int junctionIndex() { return mJunctionIndex; }
        public int extensionLength() { return mExtensionLength; }
        public int matchedBases() { return mExtensionLength - mMismatches; }
        public int mismatches() { return mMismatches; }
        public int highQualMatches() { return mHighQualMatches; }
        public boolean exceedsMaxMismatches() { return mMismatches > mPermittedMismatches; }
        public boolean exhausted() { return mExhausted; }

        public void resetMatches()
        {
            mMismatches = 0;
            mHighQualMatches = 0;
        }

        public void addMismatch() { ++mMismatches; }
        public void addHighQualMatch()
        {
            if(mCurrentIndex != mJunctionIndex)
                ++mHighQualMatches;
        }

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
            return format("%s: range(%d - %d) cigar(%s) extLen(%d) curIndex(%d) hqMatch(%d) mismatches(%d/%d) %s",
                    mRead.id(), mRead.unclippedStart(), mRead.unclippedEnd(), mRead.cigarString(),
                    mExtensionLength, mCurrentIndex, mHighQualMatches, mMismatches, mPermittedMismatches, mExhausted ? "exhausted" : "");
        }
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mBases); }
}