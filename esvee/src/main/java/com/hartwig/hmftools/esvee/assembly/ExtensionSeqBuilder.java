package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.JunctionAssembler.minReadThreshold;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findConsensusLineExtension;
import static com.hartwig.hmftools.esvee.assembly.SequenceBuilder.NEXT_BASE_CHECK_COUNT;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.common.CommonUtils.isHighBaseQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadUtils;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class ExtensionSeqBuilder
{
    private final Junction mJunction;
    private final List<ReadParseState> mReads;
    private final boolean mIsForward;
    
    private final SequenceBuilder mSequenceBuilder;

    private final int mLineExtensionLength;

    private boolean mHasLineSequence;
    private boolean mIsValid;

    // private List<RepeatInfo> mExtensionRepeats;
    // private boolean mRequiredRebuild;
    // private RepeatInfo mMaxRepeat;
    // private int mMaxRepeatCount; // includes any counts back into ref bases

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
            int readJunctionIndex = ReadUtils.getReadIndexAtJunction(read, junction, true);

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

            mReads.add(new ReadParseState(mIsForward, read, readJunctionIndex));
        }

        int baseLength = maxExtension + 1; // since the junction base itself is included (which is a ref base)
        mSequenceBuilder = new SequenceBuilder(mReads, mIsForward, baseLength);

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

        if(AssemblyConfig.AssemblyBuildDebug)
        {
            SV_LOGGER.debug("junc({}) initial extension bases sequence {}",
                    mJunction.coords(), new String(mSequenceBuilder.bases()));
        }

        validateFinalBases();
    }

    public byte[] extensionBases() { return mSequenceBuilder.bases(); }
    public int extensionLength() { return mSequenceBuilder.bases().length - 1; } // since includes the first ref base
    public byte[] baseQualities() { return mSequenceBuilder.baseQuals(); }
    public List<RepeatInfo> repeatInfo() { return mSequenceBuilder.repeats(); }
    public boolean isValid() { return mIsValid; }

    public int refBaseRepeatCount()
    {
        // TODO:
        return 0;
        // return mMaxRepeat != null ? mMaxRepeatCount - mMaxRepeat.Count : 0;
    }

    public List<SupportRead> formAssemblySupport()
    {
        List<SupportRead> supportReads = Lists.newArrayList();

        for(ReadParseState read : mReads)
        {
            if(read.mismatched())
                continue;

            if(!sufficientHighQualMatches(read))
                continue;

            supportReads.add(new SupportRead(
                    read.read(), SupportType.JUNCTION, read.startIndex(), read.matchedBases(), read.mismatches().size()));
        }

        return supportReads;
    }

    public boolean sufficientHighQualMatches(final ReadParseState read)
    {
        int refBaseRepeatBuffer = refBaseRepeatCount() > 0 ? 2 : 0;
        return read.highQualMatches() >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH + refBaseRepeatBuffer;
    }

    public List<Read> mismatchReads()
    {
        return mReads.stream().filter(x -> x.mismatched()).map(x -> x.read()).collect(Collectors.toList());
    }

    public int mismatches() { return (int)mReads.stream().filter(x -> x.mismatched()).count(); }

    public ReadParseState checkAddJunctionRead(final Read read)
    {
        int readJunctionIndex = ReadUtils.getReadIndexAtJunction(read, mJunction, true);

        if(readJunctionIndex == INVALID_INDEX)
            return null;

        // calculate how many bases beyond the junction the read extends
        // for positive orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 3
        // for negative orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
        int extensionLength = mJunction.isForward() ? read.basesLength() - readJunctionIndex - 1 : readJunctionIndex;

        ReadParseState readParseState = new ReadParseState(mIsForward, read, readJunctionIndex);

        int extensionIndex = mIsForward ? 0 : mSequenceBuilder.baseQuals().length - 1;

        int repeatIndexStart = -1;
        int repeatSkipCount = 0;

        if(mJunction.indelCoords() != null && read.indelCoords() != null
        && mJunction.indelCoords().isDelete() && read.indelCoords().isDelete())
        {
            // the read must cover the assembly INDEL entirely
            if(read.alignmentStart() < mJunction.indelCoords().PosStart && read.alignmentEnd() > mJunction.indelCoords().PosEnd
            && mJunction.indelCoords().Length != read.indelCoords().Length)
            {
                // a shorter delete means more ref bases need to be skipped
                repeatSkipCount = mJunction.indelCoords().Length - read.indelCoords().Length;
                repeatIndexStart = extensionIndex + (mJunction.isForward() ? 1 : -1); // first base after the delete
            }
        }


        final byte[] extensionBases = mSequenceBuilder.bases();
        final byte[] extensionQuals = mSequenceBuilder.baseQuals();

        while(extensionIndex >= 0 && extensionIndex < extensionBases.length)
        {
            byte consensusBase = extensionBases[extensionIndex];
            byte consensusQual = extensionQuals[extensionIndex];

            byte base = readParseState.currentBase();
            byte qual = readParseState.currentQual();
            boolean isHighQual = isHighBaseQual(qual);

            /*
            if(Nucleotides.baseIndex(base) < 0) // move to common conversion in Prep or Assembly read loading
            {
                base = DNA_N_BYTE;
                qual = INVALID_QUAL;
                isHighQual = false;
            }
            */

            if(base == consensusBase)
            {
                readParseState.addBaseMatch(isHighQual);
                readParseState.moveNext();
            }
            else
            {
                if(!isHighQual)
                {
                    // TODO: how to handle a low-qual mismatch masking a repeat or indel diff?
                }
                else
                {
                    // evaluate difference and decide how to proceed
                    assessReadMismatch(readParseState, extensionIndex, consensusBase, consensusQual);

                    if(readParseState.exhausted())
                        break;

                    if(readParseState.mismatched())
                        return readParseState;
                }
            }
        }

        if(sufficientHighQualMatches(readParseState))
        {
            // add as support
            mReads.add(readParseState);
        }

        return readParseState;
    }

    private void assessReadMismatch(final ReadParseState read, int extensionIndex, byte consensusBase, byte consensusQual)
    {
        SequenceDiffInfo seqDiffInfo = SequenceDiffInfo.UNSET;

        byte[] readNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);
        byte[] consensusNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);
        byte[] consensusCurrentAndNext = read.getNextBases(NEXT_BASE_CHECK_COUNT);

        mSequenceBuilder.populateBases(consensusCurrentAndNext, extensionIndex);
        mSequenceBuilder.populateBases(consensusNextBases, extensionIndex + (mIsForward ? 1 : -1));

        // a novel indel
        seqDiffInfo = mSequenceBuilder.assessNovelIndelDifferences(
                read, seqDiffInfo, readNextBases, consensusNextBases, consensusCurrentAndNext);

        if(seqDiffInfo == SequenceDiffInfo.UNSET)
        {
            // first an SNV match
            seqDiffInfo = SequenceDiffInfo.fromSnv(read, extensionIndex);
        }

        RepeatInfo consensusRepeat = mSequenceBuilder.findRepeat(extensionIndex);

        // TODO:
        int previousReadBases = 0;

        mSequenceBuilder.applyReadMismatches(read, seqDiffInfo, consensusRepeat, previousReadBases);
    }

    private void validateFinalBases()
    {
        int maxValidExtensionLength = 0;
        int validExtensionReadCount = 0;
        boolean hasMinLengthSoftClipRead = false;

        int reqExtensionLength = mHasLineSequence ? LINE_MIN_EXTENSION_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        int reqSecondaryExtensionLength = mHasLineSequence ? LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        for(ReadParseState read : mReads)
        {
            // CHECK, rework, make common etc
            if(read.mismatched())
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
            maxValidExtensionLength = max(read.overlapBaseCount() + 1, maxValidExtensionLength);
        }

        int minRequiredReadCount = minReadThreshold(mJunction);

        if(maxValidExtensionLength == 0 || !hasMinLengthSoftClipRead || validExtensionReadCount < minRequiredReadCount)
        {
            mIsValid = false;
        }
    }

    public String toString()
    {
        return format("junc(%s) reads(%d) baseLength(%d) lineLength(%d)",
                mJunction.coordsTyped(), mReads.size(), mSequenceBuilder.bases().length, mLineExtensionLength);
    }

    public boolean hasLineSequence() { return mHasLineSequence; }

    public String buildInformation()
    {
        // Statistics captured:
        // ReadCount - candidate read count
        // ExactMatchReads - read with no mismatches
        // MismatchReads - excluded from final support
        // Mismatch types - SNV, 1-base indels, HP repeats, and other repeats
        // Repeats
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        sj.add(String.valueOf(mReads.size()));

        int mismatchReads = (int)mReads.stream().filter(x -> x.mismatched()).count();
        int exactMatchReads = (int)mReads.stream().filter(x -> x.mismatchCount() == 0).count();
        sj.add(String.valueOf(exactMatchReads));
        sj.add(String.valueOf(mismatchReads));
        // sj.add(String.valueOf(mRequiredRebuild));

        int noRepeat = 0;
        int closeRepeat = 0;
        int notCloseRepeat = 0;
        int matchedRepeat = 0;

        /*
        if(mMaxRepeat != null)
        {
            sj.add(format("%d:%s:%d:%d", mMaxRepeat.Index, mMaxRepeat.Bases, mMaxRepeat.Count, mMaxRepeatCount - mMaxRepeat.Count));

            for(int i = 0; i < mReadRepeatCounts.length; ++i)
            {
                if(mReadRepeatCounts[i] <= 0)
                    ++noRepeat;
                else if(mReadRepeatCounts[i] == mMaxRepeat.Count)
                    ++matchedRepeat;
                else if(abs(mReadRepeatCounts[i] - mMaxRepeat.Count) <= 2)
                    ++closeRepeat;
                else
                    ++notCloseRepeat;
            }

            sj.add(format("%d:%d:%d:%d", matchedRepeat, closeRepeat, notCloseRepeat, noRepeat));
        }
        else
        {
            sj.add("0:0:0");
            sj.add("0:0:0:0");
        }
        */

        return sj.toString();
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mSequenceBuilder.bases()); }
}