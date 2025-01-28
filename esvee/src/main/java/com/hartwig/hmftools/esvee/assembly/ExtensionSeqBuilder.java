package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASE_BYTES;
import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BYTE;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.utils.Arrays.initialise;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.JunctionAssembler.minReadThreshold;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findConsensusLineExtension;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.getRepeatCount;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.getReadIndexAtReferencePosition;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class ExtensionSeqBuilder
{
    private final Junction mJunction;
    private final List<ExtReadParseState> mReads;
    private final boolean mIsForward;

    private byte[] mBases;
    private byte[] mBaseQuals;

    private final int mLineExtensionLength;

    private RepeatInfo mMaxRepeat;
    private int[] mReadRepeatCounts; // difference in number of max repeat vs the consensus
    private List<RepeatInfo> mExtensionRepeats;

    private boolean mHasLineSequence;
    private boolean mIsValid;
    private boolean mRequiredRebuild;

    private static final int READ_REPEAT_COUNT_INVALID = -1;

    public ExtensionSeqBuilder(final Junction junction, final List<Read> reads)
    {
        mJunction = junction;
        mIsForward = mJunction.isForward();
        mReads = Lists.newArrayListWithCapacity(reads.size());

        int maxExtension = 0;
        boolean hasLineReads = false;

        for(Read read : reads)
        {
            int readJunctionIndex = getReadIndexAtReferencePosition(read, junction.Position, true);

            if(readJunctionIndex == INVALID_INDEX)
                continue;

            hasLineReads |= read.hasLineTail();

            // calculate how many bases beyond the junction the read extends
            // for positive orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 3
            // for negative orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            int extensionLength = junction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

            // indel-based reads may have much shorter extensions than the typical 32-base requirement, so ensure they can extend far enough
            // to satisfy the min-high-qual match filter used post-consensus
            if(extensionLength < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                continue;

            maxExtension = max(maxExtension, extensionLength);

            mReads.add(new ExtReadParseState(read, readJunctionIndex, extensionLength, mIsForward));
        }

        int baseLength = maxExtension + 1; // since the junction base itself is included (which is a ref base)

        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];
        mExtensionRepeats = Lists.newArrayList();
        mMaxRepeat = null;
        mReadRepeatCounts = null;
        mRequiredRebuild = false;

        mIsValid = true;

        if(hasLineReads)
        {
            mLineExtensionLength = findConsensusLineExtension(reads, mJunction);
            mHasLineSequence = mLineExtensionLength >= LINE_POLY_AT_REQ;
        }
        else
        {
            mLineExtensionLength = 0;
            mHasLineSequence = false;
        }

        buildSequence(true);

        if(AssemblyConfig.AssemblyBuildDebug)
        {
            SV_LOGGER.debug("junc({}) initial extension bases sequence {}", mJunction.coords(), new String(mBases));
        }

        formConsensusSequence();

        if(AssemblyConfig.AssemblyBuildDebug)
        {
            SV_LOGGER.debug("junc({}) final extension bases sequence {}", mJunction.coords(), new String(mBases));
        }

        finaliseBases();

        findRepeats();

        assignReads();

        validateFinalBases();
    }

    public byte[] extensionBases() { return mBases; }
    public int extensionLength() { return mBases.length - 1; } // since includes the first ref base
    public byte[] baseQualities() { return mBaseQuals; }
    public List<RepeatInfo> repeatInfo() { return mExtensionRepeats; }
    public boolean isValid() { return mIsValid; }

    public List<SupportRead> formAssemblySupport()
    {
        List<SupportRead> supportReads = Lists.newArrayList();

        for(ExtReadParseState read : mReads)
        {
            if(read.exceedsMaxMismatches())
                continue;

            if(read.highQualMatches() < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
                continue;

            supportReads.add(new SupportRead(
                    read.read(), SupportType.JUNCTION, read.junctionIndex(), read.matchedBases(), read.mismatches()));
        }

        return supportReads;
    }

    public List<Read> mismatchReads()
    {
        return mReads.stream().filter(x -> x.exceedsMaxMismatches()).map(x -> x.read()).collect(Collectors.toList());
    }

    public int mismatches() { return (int)mReads.stream().filter(x -> x.exceedsMaxMismatches()).count(); }

    private void buildSequence(boolean isInitial)
    {
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        int baseCount = Nucleotides.DNA_BASES.length;

        boolean lineBasesSet = false;

        int[] readRepeatSkipCounts = new int[mReads.size()];
        int repeatIndexStart = -1; // index of the identified repeat from the junction side

        if(mMaxRepeat != null)
        {
            repeatIndexStart = mIsForward ? mMaxRepeat.Index : mMaxRepeat.postRepeatIndex() - 1;
            int repeatLength = mMaxRepeat.baseLength();
            for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
            {
                if(mReads.get(readIndex).exceedsMaxMismatches())
                    continue;

                if(mReadRepeatCounts[readIndex] != READ_REPEAT_COUNT_INVALID)
                    readRepeatSkipCounts[readIndex] = (mReadRepeatCounts[readIndex] - mMaxRepeat.Count) * repeatLength;
                else
                    readRepeatSkipCounts[readIndex] = 0;
            }
        }

        boolean checkReadRepeats = false;
        boolean readRepeatsComplete = false;
        byte[] currentReadBases = new byte[mReads.size()]; // each read's base at this index if valid and above min qual

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            byte consensusBase = 0;
            int consensusMaxQual = 0;
            int consensusQualTotal = 0;

            // per-base arrays are only used for high-qual mismatches
            int[] totalQuals = null;
            int[] maxQuals = null;

            if(mMaxRepeat != null)
            {
                if(!checkReadRepeats && !readRepeatsComplete)
                {
                    checkReadRepeats = extensionIndex == repeatIndexStart;
                }
                else if(checkReadRepeats)
                {
                    readRepeatsComplete = Arrays.stream(readRepeatSkipCounts).allMatch(x -> x == 0);
                    checkReadRepeats = !readRepeatsComplete;
                }
            }

            boolean hasActiveReads = false;

            for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
            {
                ExtReadParseState read = mReads.get(readIndex);

                currentReadBases[readIndex] = 0;

                if(read.exhausted())
                    continue;

                if(read.exceedsMaxMismatches()) // be relevant when building the final consensus
                    continue;

                hasActiveReads = true;

                byte base = read.currentBase();
                int qual = read.currentQual();

                if(aboveMinQual(qual))
                    currentReadBases[readIndex] = base;

                if(checkReadRepeats && readRepeatSkipCounts[readIndex] != 0)
                {
                    if(readRepeatSkipCounts[readIndex] < 0)
                    {
                        // this read has fewer of the current repeat than the consensus so skip using it for now
                        ++readRepeatSkipCounts[readIndex];
                        continue;
                    }
                    else if(readRepeatSkipCounts[readIndex] > 0)
                    {
                        // skip ahead of this read's extra repeats
                        read.skipBases(readRepeatSkipCounts[readIndex]);
                        readRepeatSkipCounts[readIndex] = 0;

                        if(read.exhausted())
                            break;

                        base = read.currentBase();
                        qual = read.currentQual();
                    }
                }

                if(base == DNA_N_BYTE)
                {
                    base = DNA_BASE_BYTES[0];
                    qual = 0;
                }

                // progress the read index
                read.moveNext();

                boolean hasMismatch = totalQuals != null;

                if(totalQuals == null)
                {
                    if(consensusBase == 0 || (base != consensusBase && belowMinQual(consensusMaxQual) && aboveMinQual(qual)))
                    {
                        // set first or replace with first high qual
                        consensusBase = base;
                        consensusMaxQual = qual;
                        consensusQualTotal = qual;
                    }
                    else if(base == consensusBase)
                    {
                        consensusMaxQual = max(qual, consensusMaxQual);
                        consensusQualTotal += qual;
                    }
                    else if(base != consensusBase && belowMinQual(qual))
                    {
                        // low-qual disagreement - ignore regardless of consensus qual
                    }
                    else
                    {
                        hasMismatch = true;

                        // high-qual mismatch so start tracking frequencies for each base
                        totalQuals = new int[baseCount];
                        maxQuals = new int[baseCount];

                        // back port existing counts to the per-base arrays
                        int baseIndex = Nucleotides.baseIndex(consensusBase);
                        totalQuals[baseIndex] = consensusQualTotal;
                        maxQuals[baseIndex] = consensusMaxQual;
                    }
                }

                if(hasMismatch)
                {
                    int baseIndex = Nucleotides.baseIndex(base);

                    totalQuals[baseIndex] += qual;
                    maxQuals[baseIndex] = max(maxQuals[baseIndex], qual);
                }
            }

            if(!hasActiveReads)
                break;

            if(totalQuals != null)
            {
                // take the bases with the highest qual totals
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

                if(!isInitial)
                {
                    // on a second sequence build, even after culling mismatched reads and allowing for repeat differences, some reads
                    // may still not agree - so track these differences and then stop using reads which exceed the permitted count
                    for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
                    {
                        if(currentReadBases[readIndex] > 0 && currentReadBases[readIndex] != consensusBase)
                        {
                            ExtReadParseState read = mReads.get(readIndex);
                            read.addMismatch();
                        }
                    }
                }
            }

            mBases[extensionIndex] = consensusBase;
            mBaseQuals[extensionIndex] = (byte)consensusMaxQual;

            if(mIsForward)
                ++extensionIndex;
            else
                --extensionIndex;

            if(mHasLineSequence && !lineBasesSet)
            {
                setLineExtensionBases(extensionIndex);
                lineBasesSet = true;

                if(mIsForward)
                    extensionIndex += mLineExtensionLength;
                else
                    extensionIndex -= mLineExtensionLength;
            }
        }
    }

    private byte lineBase() { return mJunction.isForward() ? LINE_BASE_T : LINE_BASE_A; }

    private void setLineExtensionBases(int extensionIndex)
    {
        // build out line bases if identified
        byte lineBase = lineBase();

        int remainingBases = mLineExtensionLength;

        while(extensionIndex >= 0 && extensionIndex < mBases.length && remainingBases > 0)
        {
            mBases[extensionIndex] = lineBase;
            mBaseQuals[extensionIndex] = (byte)LOW_BASE_QUAL_THRESHOLD;

            extensionIndex += mIsForward ? 1 : -1;

            --remainingBases;
        }

        // move each read to the end of its poly A/T sequence
        moveReadsPastLineExtension();
    }

    private void moveReadsPastLineExtension()
    {
        byte lineBase = lineBase();
        mReads.forEach(x -> x.movePastLineExtension(lineBase, false));
    }

    public void findConsensusRepeats()
    {
        if(mHasLineSequence)
            return;

        // check the existing sequence for a valid repeat
        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);

        if(repeats.isEmpty())
            return;

        mReadRepeatCounts = new int[mReads.size()];

        Collections.sort(repeats, Comparator.comparingInt(x -> -x.Count));
        RepeatInfo maxRepeat = repeats.get(0);

        // look for a common repeat with variation from all reads, indexed from the junction
        int repeatIndexStart = mIsForward ? maxRepeat.Index : maxRepeat.postRepeatIndex() - 1;
        int repeatJunctionOffset = mIsForward ? repeatIndexStart : mBases.length - 1 - repeatIndexStart;

        Map<Integer,Integer> skipFrequencies = Maps.newHashMap();
        List<Integer> readRepeatCounts = Lists.newArrayList();

        for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
        {
            ExtReadParseState read = mReads.get(readIndex);

            int readJunctionIndex = read.junctionIndex();
            int readRepeatIndexStart = readJunctionIndex + (mIsForward ? repeatJunctionOffset : -repeatJunctionOffset);

            int readRepeatCount = getRepeatCount(read.read(), maxRepeat, readRepeatIndexStart, mIsForward);

            if(readRepeatCount >= 0) // -1 means the repeat section is outside the bounds of this read, so ignore it for frequencies
            {
                mReadRepeatCounts[readIndex] = readRepeatCount;
                Integer freq = skipFrequencies.get(readRepeatCount);
                skipFrequencies.put(readRepeatCount, freq != null ? freq + 1 : 1);
                readRepeatCounts.add(readRepeatCount);
            }
            else
            {
                mReadRepeatCounts[readIndex] = READ_REPEAT_COUNT_INVALID;
            }
        }

        if(skipFrequencies.size() == 1)
            return; // no disagreement, mismatches must be explained by other read sequence differences

        // first check the median repeat count and if zero then assume there is no dominant repeats
        int medianRepeatCount = Doubles.medianInteger(readRepeatCounts);

        if(medianRepeatCount == 0)
            return;

        int consensusRepeatCount = 0;
        int consensusRepeatFreq = 0;

        for(Map.Entry<Integer,Integer> entry : skipFrequencies.entrySet())
        {
            if(entry.getKey() == 0) // ignore zero has the most frequent since the median is above zero
                continue;

            if(entry.getValue() > consensusRepeatFreq)
            {
                consensusRepeatFreq = entry.getValue();
                consensusRepeatCount = entry.getKey();
            }
        }

        // recompute the median from the non-zero repeat counts
        readRepeatCounts = readRepeatCounts.stream().filter(x -> x.intValue() > 0).collect(Collectors.toList());
        medianRepeatCount = Doubles.medianInteger(readRepeatCounts);

        // use median if there is no dominant repeat count
        if(consensusRepeatFreq < readRepeatCounts.size() / 2)
            consensusRepeatCount = medianRepeatCount;

        repeatIndexStart = maxRepeat.Index;

        // adjust the repeat start index on the reverse strand to the consensus count
        if(!mIsForward && consensusRepeatCount != maxRepeat.Count)
            repeatIndexStart -= (consensusRepeatCount - maxRepeat.Count) * maxRepeat.baseLength();

        mMaxRepeat = new RepeatInfo(repeatIndexStart, maxRepeat.Bases, consensusRepeatCount);
    }

    private void formConsensusSequence()
    {
        // test reads against the initial sequence, filtering out those with too many mismatches, and forming a consensus from the rest
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        mReads.forEach(x -> x.resetIndex());

        if(mHasLineSequence)
        {
            extensionIndex += mIsForward ? (mLineExtensionLength + 1) : -(mLineExtensionLength + 1);
            moveReadsPastLineExtension();
        }

        int[] readMismatchExceededIndex = new int[mReads.size()];
        initialise(readMismatchExceededIndex, -1);

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
            {
                ExtReadParseState read = mReads.get(readIndex);

                if(read.exhausted())
                    continue;

                byte base = read.currentBase();
                byte qual = read.currentQual();

                if(aboveMinQual(qual) && aboveMinQual(mBaseQuals[extensionIndex]))
                {
                    if(base != mBases[extensionIndex])
                    {
                        read.addMismatch();

                        // record the first base where the read exceeds its permitted mismatch count
                        if(read.exceedsMaxMismatches() && readMismatchExceededIndex[readIndex] < 0)
                        {
                            readMismatchExceededIndex[readIndex] = extensionIndex;
                        }
                    }
                    else
                    {
                        read.addHighQualMatch();
                    }
                }

                read.moveNext();
            }

            extensionIndex += mIsForward ? 1 : -1;
        }

        int mismatchReadCount = (int)mReads.stream().filter(x -> x.exceedsMaxMismatches()).count();

        // if no reads have been excluded then keep the initial sequence
        if(mismatchReadCount == 0)
            return;

        mRequiredRebuild = true;

        // otherwise look for a variable repeat as an explanation for the mismatches
        findConsensusRepeats();

        // and if one is found, reset any read whose mismatches might be explained by the repeat
        if(mMaxRepeat != null)
        {
            for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
            {
                ExtReadParseState read = mReads.get(readIndex);

                if(!read.exceedsMaxMismatches() || mReadRepeatCounts[readIndex] == mMaxRepeat.Count)
                {
                    read.resetMatches();
                    continue;
                }

                // read must differ vs the consensus repeat
                if(mReadRepeatCounts[readIndex] == READ_REPEAT_COUNT_INVALID)
                    continue;

                // the read must have any at least one repeat to be considered to differ from jitter
                if(mReadRepeatCounts[readIndex] == 0 && mMaxRepeat.Count > 1)
                    continue;

                // if a read has repeats and its mismatch occurs within the range of the first repeat, consider it mismatched from jittter
                int mismatchExceededIndex = readMismatchExceededIndex[readIndex];

                int firstRepeatEndIndex = mIsForward ?
                        mMaxRepeat.Index + mMaxRepeat.baseLength() : mMaxRepeat.postRepeatIndex() - 1 - mMaxRepeat.baseLength();

                if(mIsForward && mismatchExceededIndex >= firstRepeatEndIndex)
                    read.resetMatches();
                else if(!mIsForward && mismatchExceededIndex <= firstRepeatEndIndex)
                    read.resetMatches();
            }
        }

        // otherwise filter out the mismatch reads and build the sequence again
        int baseLength = mBases.length;
        mBases = new byte[baseLength];
        mBaseQuals = new byte[baseLength];

        mReads.forEach(x -> x.resetIndex());

        buildSequence(false);
    }

    private void finaliseBases()
    {
        // trim extension bases if required
        int validBaseCount = 0;
        int index = mIsForward ? 0 : mBases.length - 1;

        while(index >= 0 && index < mBases.length)
        {
            if(mBases[index] == 0)
                break;

            ++validBaseCount;
            index += mIsForward ? 1 : -1;
        }

        if(validBaseCount == mBases.length)
            return;

        int reduction = mBases.length - validBaseCount;
        int startIndex = mIsForward ? 0 : reduction;
        int endIndex = mIsForward ? validBaseCount - 1 : mBases.length - 1;
        mBases = subsetArray(mBases, startIndex, endIndex);
        mBaseQuals = subsetArray(mBaseQuals, startIndex, endIndex);

        if(!mIsForward && mMaxRepeat != null)
            mMaxRepeat = new RepeatInfo(mMaxRepeat.Index - reduction, mMaxRepeat.Bases, mMaxRepeat.Count);
    }

    private void findRepeats()
    {
        mExtensionRepeats.clear();

        List<RepeatInfo> repeats = RepeatInfo.findRepeats(mBases);

        if(repeats != null)
            mExtensionRepeats.addAll(repeats);
    }

    private void assignReads()
    {
        // check all reads against the final sequence with awareness of LINE and any consensus repeat
        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        if(mHasLineSequence)
            extensionIndex += mIsForward ? (mLineExtensionLength + 1) : -(mLineExtensionLength + 1);

        int repeatIndexStart = -1;
        int repeatLength = 0;

        if(mMaxRepeat != null)
        {
            repeatIndexStart = mIsForward ? mMaxRepeat.Index : mMaxRepeat.postRepeatIndex() - 1;
            repeatLength = mMaxRepeat.baseLength();
        }

        for(int readIndex = 0; readIndex < mReads.size(); ++readIndex)
        {
            ExtReadParseState read = mReads.get(readIndex);

            read.resetIndex();
            read.resetMatches();

            int readRepeatSkipCount = mMaxRepeat != null && mReadRepeatCounts[readIndex] != READ_REPEAT_COUNT_INVALID ?
                    (mReadRepeatCounts[readIndex] - mMaxRepeat.Count) * repeatLength : 0;

            checkReadMismatches(read, extensionIndex, repeatIndexStart, readRepeatSkipCount);
        }
    }

    public ExtReadParseState checkAddJunctionRead(final Read read)
    {
        int readJunctionIndex = getReadIndexAtReferencePosition(read, mJunction.Position, true);

        if(readJunctionIndex == INVALID_INDEX)
            return null;

        // calculate how many bases beyond the junction the read extends
        // for positive orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 3
        // for negative orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
        int extensionLength = mJunction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

        ExtReadParseState readParseState = new ExtReadParseState(read, readJunctionIndex, extensionLength, mIsForward);

        int extensionIndex = mIsForward ? 0 : mBases.length - 1;

        if(mHasLineSequence)
            extensionIndex += mIsForward ? (mLineExtensionLength + 1) : -(mLineExtensionLength + 1);

        int repeatIndexStart = -1;
        int repeatSkipCount = 0;

        if(mMaxRepeat != null)
        {
            repeatIndexStart = mIsForward ? mMaxRepeat.Index : mMaxRepeat.postRepeatIndex() - 1;

            int repeatJunctionOffset = mIsForward ? repeatIndexStart : mBases.length - 1 - repeatIndexStart;
            int readRepeatIndexStart = readJunctionIndex + (mIsForward ? repeatJunctionOffset : -repeatJunctionOffset);
            int readRepeatCount = getRepeatCount(read, mMaxRepeat, readRepeatIndexStart, mIsForward);

            repeatSkipCount = (readRepeatCount - mMaxRepeat.Count) * mMaxRepeat.baseLength();
        }

        checkReadMismatches(readParseState, extensionIndex, repeatIndexStart, repeatSkipCount);

        return readParseState;
    }

    private void checkReadMismatches(final ExtReadParseState read, int extensionIndex, int repeatIndexStart, final int repeatSkipCount)
    {
        if(mHasLineSequence)
            read.movePastLineExtension(lineBase(), true);

        int remainingRepeatSkipCount = repeatSkipCount;
        boolean checkReadRepeats = false;

        while(extensionIndex >= 0 && extensionIndex < mBases.length)
        {
            if(read.exhausted() || read.exceedsMaxMismatches())
                break;

            byte base = read.currentBase();
            byte qual = read.currentQual();

            if(extensionIndex == repeatIndexStart)
                checkReadRepeats = true;

            if(checkReadRepeats && remainingRepeatSkipCount != 0)
            {
                if(remainingRepeatSkipCount < 0)
                {
                    // this read has fewer of the current repeat than the consensus so skip using it for now
                    ++remainingRepeatSkipCount;
                    extensionIndex += mIsForward ? 1 : -1;
                    continue;
                }
                else if(remainingRepeatSkipCount > 0)
                {
                    // skip ahead of this read's extra repeats
                    read.skipBases(remainingRepeatSkipCount);
                    remainingRepeatSkipCount = 0;

                    if(read.exhausted())
                        break;

                    base = read.currentBase();
                    qual = read.currentQual();
                }
            }

            if(aboveMinQual(qual) && aboveMinQual(mBaseQuals[extensionIndex]))
            {
                if(base != mBases[extensionIndex])
                    read.addMismatch();
                else
                    read.addHighQualMatch();
            }

            read.moveNext();

            extensionIndex += mIsForward ? 1 : -1;
        }
    }

    private void validateFinalBases()
    {
        int maxValidExtensionLength = 0;
        int validExtensionReadCount = 0;
        boolean hasMinLengthSoftClipRead = false;

        int reqExtensionLength = mHasLineSequence ? LINE_MIN_EXTENSION_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        int reqSecondaryExtensionLength = mHasLineSequence ? LINE_MIN_EXTENSION_LENGTH / 2 : ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        for(ExtReadParseState read : mReads)
        {
            if(read.exceedsMaxMismatches())
                continue;

            int scLength;

            if(mJunction.DiscordantOnly)
            {
                scLength = mJunction.isForward() ?
                        read.read().unclippedEnd() - mJunction.Position : mJunction.Position - read.read().unclippedStart();
            }
            else
            {
                scLength = mJunction.isForward() ? read.read().rightClipLength() : read.read().leftClipLength();
            }

            if(scLength >= reqSecondaryExtensionLength)
            {
                hasMinLengthSoftClipRead |= scLength >= reqExtensionLength;
                ++validExtensionReadCount;
            }
            else if(read.read().indelCoords() != null && read.read().indelCoords().Length >= MIN_INDEL_LENGTH)
            {
                // check for extensions comprised only of short indels, even if they started with some soft-clipped reads
                ++validExtensionReadCount;
                hasMinLengthSoftClipRead = true;
            }
            else
            {
                continue;
            }

            // note the read's extension length does not include the junction base itself, hence the +1
            maxValidExtensionLength = max(read.extensionLength() + 1, maxValidExtensionLength);
        }

        int minRequiredReadCount = minReadThreshold(mJunction);

        if(maxValidExtensionLength == 0 || !hasMinLengthSoftClipRead || validExtensionReadCount < minRequiredReadCount)
        {
            mIsValid = false;
        }
    }

    public String toString()
    {
        return format("junc(%s) reads(%d) baseLength(%d) lineLength(%d) maxRepeat(%s)",
                mJunction.coordsTyped(), mReads.size(), mBases.length, mLineExtensionLength,
                mMaxRepeat != null ? mMaxRepeat : "none");
    }

    public boolean hasLineSequence() { return mHasLineSequence; }

    public String buildInformation()
    {
        // ReadCount;MismatchReadCount;Rebuild;MaxRepeatInfo;ReadRepeatInfo;
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        sj.add(String.valueOf(mReads.size()));

        int mismatchReads = (int)mReads.stream().filter(x -> x.exceedsMaxMismatches()).count();
        sj.add(String.valueOf(mismatchReads));
        sj.add(String.valueOf(mRequiredRebuild));

        int noRepeat = 0;
        int closeRepeat = 0;
        int notCloseRepeat = 0;
        int matchedRepeat = 0;

        if(mMaxRepeat != null)
        {
            sj.add(format("%d:%s:%d", mMaxRepeat.Index, mMaxRepeat.Bases, mMaxRepeat.Count));

            for(int i = 0; i < mReadRepeatCounts.length; ++i)
            {
                if(mReadRepeatCounts[i] <= 0)
                    ++noRepeat;
                else if(mReadRepeatCounts[i] == mMaxRepeat.Count)
                    ++matchedRepeat;
                else if(abs(mReadRepeatCounts[i] - mMaxRepeat.Count) <= 2)
                    ++closeRepeat;
                else if(mReadRepeatCounts[i] == mMaxRepeat.Count)
                    ++notCloseRepeat;
            }

            sj.add(format("%d:%d:%d:%d", matchedRepeat, closeRepeat, notCloseRepeat, noRepeat));
        }
        else
        {
            sj.add("0:0:0");
            sj.add("0:0:0:0");
        }
        return sj.toString();
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mBases); }
}