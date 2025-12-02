package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.JunctionAssembler.minReadThreshold;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findConsensusLineExtension;
import static com.hartwig.hmftools.esvee.assembly.SequenceBuilder.NEXT_BASE_CHECK_COUNT;
import static com.hartwig.hmftools.esvee.assembly.SequenceBuilder.getRepeatCount;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.BASE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.DELETE;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.INSERT;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.REPEAT;
import static com.hartwig.hmftools.esvee.assembly.SequenceDiffType.UNSET;
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
    private final boolean mBuildForwards;
    
    private final SequenceBuilder mSequenceBuilder;

    private boolean mHasLineSequence;
    private boolean mIsValid;

    private int mRefRepeatCount; // counts of a repeat which extends back into the ref bases

    public ExtensionSeqBuilder(final Junction junction, final List<Read> reads)
    {
        mJunction = junction;
        mBuildForwards = mJunction.isForward();
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

            mReads.add(new ReadParseState(mBuildForwards, read, readJunctionIndex));
        }

        int baseLength = maxExtension + 1; // since the junction base itself is included (which is a ref base)
        mSequenceBuilder = new SequenceBuilder(mReads, mBuildForwards, baseLength);

        mIsValid = true;
        mRefRepeatCount = 0;

        if(hasLineReads)
        {
            int lineExtensionLength = findConsensusLineExtension(reads, mJunction);
            mHasLineSequence = lineExtensionLength >= LINE_POLY_AT_REQ;
        }
        else
        {
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
    public List<RepeatInfo> repeats() { return mSequenceBuilder.repeats(); }
    public boolean isValid() { return mIsValid; }
    public int refBaseRepeatCount() { return mRefRepeatCount; }

    public List<SupportRead> formAssemblySupport()
    {
        List<SupportRead> supportReads = Lists.newArrayList();

        for(ReadParseState read : mReads)
        {
            if(read.mismatched())
                continue;

            if(!sufficientQualMatches(read))
                continue;

            SupportRead supportRead = new SupportRead(
                    read.read(), SupportType.JUNCTION, read.startIndex(), read.matchedBases(), read.mismatchCount(true));
            supportRead.setMismatchInfo(read.mismatchInfo());

            supportReads.add(supportRead);
        }

        return supportReads;
    }

    public boolean sufficientQualMatches(final ReadParseState read)
    {
        int refBaseRepeatBuffer = refBaseRepeatCount() > 0 ? 2 : 0;
        int highQualExtensionMatches = read.highQualMatches() - 1; // to remove the initial junction ref bases
        return highQualExtensionMatches >= ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH + refBaseRepeatBuffer;
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

        ReadParseState readParseState = new ReadParseState(mBuildForwards, read, readJunctionIndex);

        int extensionIndex = mBuildForwards ? 0 : mSequenceBuilder.baseQuals().length - 1;

        int repeatIndexStart = -1;
        int repeatSkipCount = 0;

        // TODO: check how this needs to work for reads with shorter DELs
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

        while(extensionIndex >= 0 && extensionIndex < extensionBases.length)
        {
            byte consensusBase = extensionBases[extensionIndex];

            byte base = readParseState.currentBase();
            byte qual = readParseState.currentQual();
            boolean isHighQual = isHighBaseQual(qual);

            int extBaseMove = 1;

            if(base == consensusBase)
            {
                readParseState.addBaseMatch(isHighQual);
                readParseState.moveNext();
            }
            else
            {
                // evaluate difference and decide how to proceed in the read and extension bases
                extBaseMove = assessReadMismatch(readParseState, extensionIndex);

                if(readParseState.mismatched())
                    return readParseState;
            }

            if(readParseState.exhausted())
                break;

            extensionIndex += mBuildForwards ? extBaseMove : -extBaseMove;
        }

        if(sufficientQualMatches(readParseState))
        {
            // add as support
            mReads.add(readParseState);
        }

        return readParseState;
    }

    private int assessReadMismatch(final ReadParseState read, int extensionIndex)
    {
        SequenceDiffInfo seqDiffInfo = SequenceDiffInfo.UNSET;

        byte[] readNextBases = read.getNextBases(NEXT_BASE_CHECK_COUNT);
        byte[] consensusNextBases = new byte[NEXT_BASE_CHECK_COUNT];
        byte[] consensusCurrentAndNext = new byte[NEXT_BASE_CHECK_COUNT];

        mSequenceBuilder.populateBases(consensusCurrentAndNext, extensionIndex);
        mSequenceBuilder.populateBases(consensusNextBases, extensionIndex + (mBuildForwards ? 1 : -1));

        // follow the same ordering for mismatch types as in the sequence building routine:
        // SNVs
        // homopolymer repeat diffs
        // indels (limited to a single base for now)
        // longer repeat diffs

        RepeatInfo consensusRepeat = null;
        int previousRepeatLength = 0;

        if(readNextBases == null || consensusNextBases == null || basesMatch(readNextBases, consensusNextBases))
        {
            // this will handle the low-qual mismatch
            seqDiffInfo = SequenceDiffInfo.fromSnv(read, extensionIndex);
        }
        else
        {
            // first check for repeats

            // for the read to be a possible repeat adjustment (ie expansion or subtraction), it must have a repeat starting on the same
            // base as the consensus repeat - so first check if such a consensus repeat exists that would reach this base

            int priorExtIndex = extensionIndex + (mBuildForwards ? -1 : 1);
            consensusRepeat = mSequenceBuilder.findRepeat(priorExtIndex);

            if(consensusRepeat != null)
            {
                int prevExtBaseLength = abs(consensusRepeat.Index - extensionIndex); // was +1 but seems incorrect
                // previousRepeatLength = abs(consensusRepeat.Index - extensionIndex) + 1;

                int repeatBaseLength = consensusRepeat.repeatLength();
                int previousExtRepeatCount = prevExtBaseLength / repeatBaseLength;
                int repeatExtRemainder = prevExtBaseLength % repeatBaseLength;

                int readRepeatCount = 0;

                if(previousExtRepeatCount > 0 && repeatExtRemainder == 0)
                {
                    // get further repeat counts in the read from this point, assuming has matched until this point
                    readRepeatCount = getRepeatCount(read, consensusRepeat.Bases, previousExtRepeatCount, mBuildForwards);
                }

                if(readRepeatCount > 0)
                {
                    previousRepeatLength = prevExtBaseLength;

                    int readRepeatStart = SequenceDiffInfo.repeatIndex(
                            read.readIndex(), readRepeatCount, consensusRepeat, previousRepeatLength, mBuildForwards, true);

                    int readRepeatEnd = SequenceDiffInfo.repeatIndex(
                            read.readIndex(), readRepeatCount, consensusRepeat, previousRepeatLength, mBuildForwards, false);

                    int readRepeatIndexBegin = mBuildForwards ? readRepeatStart : readRepeatEnd;

                    BaseQualType baseQualType = read.rangeMinQualType(readRepeatStart, readRepeatEnd);

                    seqDiffInfo = new SequenceDiffInfo(
                            read.readIndex(), extensionIndex, consensusRepeat.Bases, REPEAT, baseQualType,
                            readRepeatCount, readRepeatIndexBegin, 0);
                }
            }

            if(seqDiffInfo.Type == UNSET)
            {
                seqDiffInfo = mSequenceBuilder.assessNovelIndelDifferences(
                        read, seqDiffInfo, readNextBases, consensusNextBases, consensusCurrentAndNext);
            }
       }

        if(seqDiffInfo == SequenceDiffInfo.UNSET)
        {
            // first an SNV match
            seqDiffInfo = SequenceDiffInfo.fromSnv(read, extensionIndex);
        }

        mSequenceBuilder.applyReadMismatches(read, seqDiffInfo, consensusRepeat, previousRepeatLength);

        int extBaseMove = 1;

        // NOTE: indels don't need to change the extension index move since the next bases (including current) have been checked for a match
        if(seqDiffInfo.Type == REPEAT)
        {
            // if the consensus repeat is longer than the read's repeat, then the extension index will move further ahead
            // otherwise it will stay as is, waiting to compare the next base again
            if(consensusRepeat.Count > seqDiffInfo.RepeatCount)
                extBaseMove = (consensusRepeat.Count - seqDiffInfo.RepeatCount) * consensusRepeat.repeatLength();
            else if(seqDiffInfo.RepeatCount > consensusRepeat.Count)
                extBaseMove = 0;
        }

        return extBaseMove;
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
            return;
        }

        checkRepeatInRefBases();
    }

    private void checkRepeatInRefBases()
    {
        if(mSequenceBuilder.repeats().isEmpty())
            return;

        RepeatInfo firstRepeat = mSequenceBuilder.repeats().get(0);

        // check if the repeat starts at or close enough to the start of the junction to be considered
        int repeatJunctionDistance = mBuildForwards ? firstRepeat.Index : mSequenceBuilder.baseLength() - firstRepeat.Index - 1;

        if(repeatJunctionDistance > firstRepeat.repeatLength())
            return;

        ReadParseState sampleRead = mReads.stream().filter(x -> !x.mismatched()).findFirst().orElse(null);

        if(sampleRead == null)
            return;

        int junctionIndex = sampleRead.startIndex();

        int repeatIndexStart = junctionIndex + (mBuildForwards ? repeatJunctionDistance - 1 : -repeatJunctionDistance + 1);
        int previousRepeatCount = RepeatInfo.getRepeatCount(sampleRead.read().getBases(), firstRepeat.Bases, repeatIndexStart, !mBuildForwards);

        if(previousRepeatCount > 0)
            mRefRepeatCount = previousRepeatCount;
    }

    public String toString()
    {
        return format("junc(%s) reads(%d) baseLength(%d)", mJunction.coordsTyped(), mReads.size(), mSequenceBuilder.bases().length);
    }

    public boolean hasLineSequence() { return mHasLineSequence; }

    public String buildInformation()
    {
        // Statistics captured:

        // RC: ReadCount - candidate read count
        // EX: ExactMatchReads - read with no mismatches
        // LQ: low qual mismatches
        // MM: MismatchReads - excluded from final support
        // SNV: SNV counts excluding low-qual
        // HPR: homopolymer repeats
        // OR: other repeats
        // ID: indels
        StringJoiner sj = new StringJoiner(ITEM_DELIM);

        int exactMatch = 0;
        int mismatchedReads = 0;
        int lowQualMismatches = 0;
        int snvs = 0;
        int indels = 0;
        int homopolymers = 0;
        int otherRepeats = 0;

        for(ReadParseState read : mReads)
        {
            if(!sufficientQualMatches(read))
                continue;

            if(read.mismatchCount(false) == 0)
            {
                ++exactMatch;
            }
            else if(read.mismatched())
            {
                ++mismatchedReads;
            }
            else
            {
                for(SequenceDiffInfo mismatch : read.mismatches())
                {
                    if(mismatch.MismatchPenalty > 0)
                    {
                        if(mismatch.Type == BASE)
                        {
                            ++snvs;
                        }
                        else if(mismatch.Type == REPEAT)
                        {
                            if(mismatch.Bases.length() == 1)
                                ++homopolymers;
                            else
                                ++otherRepeats;
                        }
                        else if(mismatch.Type == INSERT || (mismatch.Type == DELETE))
                        {
                            ++indels;
                        }
                    }
                    else
                    {
                        ++lowQualMismatches;
                    }
                }
            }
        }

        sj.add(format("RC=%d", mReads.size()));
        sj.add(format("EM=%d", exactMatch));
        sj.add(format("LQ=%d", lowQualMismatches));
        sj.add(format("MM=%d", mismatchedReads));
        sj.add(format("SNV=%d", snvs));
        sj.add(format("HP=%d", homopolymers));
        sj.add(format("OR=%d", otherRepeats));
        sj.add(format("ID=%d", indels));

        for(RepeatInfo repeat : mSequenceBuilder.repeats())
        {
            sj.add(format("REP:%dx%s", repeat.Count, repeat.Bases));
        }

        return sj.toString();
    }

    @VisibleForTesting
    public String junctionSequence() { return new String(mSequenceBuilder.bases()); }
}