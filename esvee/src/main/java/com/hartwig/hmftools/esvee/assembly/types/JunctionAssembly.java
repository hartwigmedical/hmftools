package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.alignment.AlignmentOutcome.NO_SET;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.readQualFromJunction;
import static com.hartwig.hmftools.esvee.assembly.ExtensionSeqBuilder.calcReadSequenceMismatches;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.convertedIndelCrossesJunction;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findInsertedBases;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.LOCAL_REF;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.SHORT_INDEL;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.REMOTE_REGION;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.findUnsetBases;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.calcTrimmedRefBaseLength;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;

import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.alignment.AlignmentOutcome;
import com.hartwig.hmftools.esvee.assembly.filters.FilterType;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class JunctionAssembly
{
    private int mAssemblyId;
    private final Junction mJunction;

    private int mJunctionSequenceIndex; // position of the junction in the read bases

    // aligned position on the extension side is an inferred position only
    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private final IndelCoords mIndelCoords;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<SupportRead> mSupport;
    private final List<SupportRead> mCandidateSupport;
    private final List<RefSideSoftClip> mRefSideSoftClips;

    private final List<RepeatInfo> mRepeatInfo;
    private String mRefBasesRepeatedTrimmed;
    private int mRefBaseTrimLength;

    private final List<RemoteRegion> mRemoteRegions;

    private PhaseGroup mPhaseGroup;
    private AssemblyOutcome mOutcome;
    private AlignmentOutcome mAlignmentOutcome;

    // info only
    private final String mInitialReadId;
    private int mMergedAssemblies;
    private int mMismatchReadCount;

    private final Set<FilterType> mFilters;
    private final AssemblyStats mStats;

    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] baseQualities, final List<SupportRead> assemblySupport,
            final List<RepeatInfo> repeatInfo)
    {
        mJunction = junction;

        // set initial bounds from the assembly support
        int minAlignedPosition = 0;
        int maxAlignedPosition = 0;
        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        IndelCoords indelCoords = null;
        String indelInsertedBases = "";

        // FIXME: consider moving this out of the constructor and set prior - removes logic and possible dependence on Read being cached
        for(SupportRead support : assemblySupport)
        {
            Read read = support.cachedRead();

            maxAlignedPosition = max(maxAlignedPosition, read.alignmentEnd());
            minAlignedPosition = minAlignedPosition == 0 ? read.alignmentStart() : min(minAlignedPosition, read.alignmentStart());

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }

            if(read.indelCoords() != null && indelCoords == null)
            {
                indelCoords = read.indelCoords();

                if(indelCoords.isInsert())
                    indelCoords.setInsertedBases(findInsertedBases(read));
            }
        }

        mIndelCoords = indelCoords;

        mInitialReadId = maxJunctionBaseQualRead != null ? maxJunctionBaseQualRead.id() :
                (!assemblySupport.isEmpty() ? assemblySupport.get(0).id() : "null");

        mBases = bases;
        mBaseQuals = baseQualities;

        mJunctionSequenceIndex = junction.isForward() ? 0 : mBases.length - 1;
        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        mSupport = Lists.newArrayList(assemblySupport);
        mCandidateSupport = Lists.newArrayList();
        mRepeatInfo = repeatInfo;
        mRefSideSoftClips = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mMergedAssemblies = 0;
        mPhaseGroup = null;
        mOutcome = UNSET;
        mAlignmentOutcome = NO_SET;
        mFilters = Sets.newHashSet();
        mMismatchReadCount = 0;
        mStats = new AssemblyStats();
    }

    public void setId(int id) { mAssemblyId = id; }
    public int id() { return mAssemblyId; }

    public Junction junction() { return mJunction; }
    public boolean isForwardJunction() { return mJunction.isForward(); }

    public boolean indel() { return mJunction.IndelBased; }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionSequenceIndex; };
    public void setJunctionIndex(int index) { mJunctionSequenceIndex = index; };

    // eg 21 bases, junction index at 10 (so 0-9 = 10 before, 11-20 = 10 after), note: doesn't count the junction base
    public int lowerDistanceFromJunction() { return mJunctionSequenceIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionSequenceIndex - 1; };

    public int refBaseLength() { return (mJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction()) + 1; }

    // base count beyond the junction
    public int extensionLength() { return mJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    // purely for informational purposes at this stage - since past the junction is just the soft-clip length, and on the
    // reference side doesn't take into account any INDELs
    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int refBasePosition() { return mJunction.isForward() ? mMinAlignedPosition : mMaxAlignedPosition; }
    public int baseLength() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public String initialReadId() { return mInitialReadId; }
    public IndelCoords indelCoords() { return mIndelCoords; }

    public List<SupportRead> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public AssemblyStats stats() { return mStats; }

    public boolean checkAddJunctionRead(final Read read, int permittedMismatches)
    {
        int mismatchCount = 0;

        int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);

        if(readJunctionIndex == INVALID_INDEX)
            return false;

        int[] readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, false);
        int assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], false);

        if(assemblyIndex < 0)
            return false;

        int highQualMatchCount = 0;

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length)
                break;

            byte qual = read.getBaseQuality()[i];

            if(basesMatch(read.getBases()[i], mBases[assemblyIndex], qual, mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                if(qual >= LOW_BASE_QUAL_THRESHOLD)
                    ++highQualMatchCount;
            }
            else
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    break;
            }
        }

        mismatchCount = calcReadSequenceMismatches(
                mJunction.isForward(), mBases, mBaseQuals, mRepeatInfo, read, readJunctionIndex, permittedMismatches);

        if(mismatchCount > permittedMismatches)
            return false;

        addSupport(read, JUNCTION, readJunctionIndex, highQualMatchCount, mismatchCount);
        return true;
    }

    public int mismatchReadCount() { return mMismatchReadCount; }
    public void addMismatchReadCount(int count) { mMismatchReadCount += count; }

    public void extendJunctionReadSupport(final Read read, final SupportRead existingSupport)
    {
        addRead(read, existingSupport.type(), existingSupport);
    }

    private void addRead(final Read read, final SupportType type, @Nullable final SupportRead existingSupport)
    {
        int mismatchCount = 0;
        int highQualMatchCount = 0;
        int readJunctionIndex;
        int[] readIndexRange;
        int assemblyIndex;

        boolean byReference = existingSupport != null;

        if(type == JUNCTION || type == INDEL)
        {
            if(existingSupport == null)
            {
                readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);

                if(readJunctionIndex == INVALID_INDEX)
                    return;
            }
            else
            {
                readJunctionIndex = existingSupport.junctionReadIndex();
            }

            readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, byReference);

            assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], byReference);
        }
        else
        {
            readIndexRange = new int[] { 0, read.getBases().length - 1 }; // take the whole read for ref-side reads
            readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, false);

            if(mJunction.isForward())
            {
                assemblyIndex = read.unclippedStart() - mMinAlignedPosition;
            }
            else
            {
                assemblyIndex = mJunctionSequenceIndex + read.unclippedStart() - mJunction.Position;
            }
        }

        if(readIndexRange[SE_START] < 0)
            return;

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length || i >= read.getBases().length)
            {
                // can possibly happen with INDELs in the CIGAR or not now that aligned positions are used to establish read coordinates?
                // SV_LOGGER.debug("i({}) readCoords({}) assembly({}) read({})", i, readCoords, toString(), read.toString());
                break;
            }

            byte base = read.getBases()[i];
            byte qual = read.getBaseQuality()[i];

            if(mBases[assemblyIndex] == 0)
            {
                mBases[assemblyIndex] = base;
                mBaseQuals[assemblyIndex] = qual;

                if(qual >= LOW_BASE_QUAL_THRESHOLD)
                    ++highQualMatchCount;
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    if(qual >= LOW_BASE_QUAL_THRESHOLD)
                        ++highQualMatchCount;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                }
                else
                {
                    ++mismatchCount;
                }
            }
        }

        if(existingSupport == null)
        {
            addSupport(read, type, readJunctionIndex, highQualMatchCount, mismatchCount);
        }
        else
        {
            // not set since not used for now
            //int[] existingReadRange = existingSupport.readIndexRange();
            // existingSupport.setReadIndexRange(
            //        min(existingReadRange[SE_START], readIndexRange[SE_START]), max(existingReadRange[SE_END], readIndexRange[SE_END]));

            existingSupport.setReferenceMismatches(mismatchCount);
        }
    }

    private void addSupport(final Read read, final SupportType type, final int readJunctionIndex, final int matches, final int mismatches)
    {
        boolean isIndelCrossingJunction = convertedIndelCrossesJunction(mJunction, read);
        SupportType adjustedType = type == JUNCTION && isIndelCrossingJunction ? INDEL : type;
        SupportRead support = new SupportRead(read, adjustedType, readJunctionIndex, matches, mismatches);

        mSupport.add(support);

        mStats.addRead(support, mJunction, read);

        if(mOutcome != SHORT_INDEL && !indel() && adjustedType == INDEL)
        {
            setOutcome(SHORT_INDEL);
        }
    }

    public void extendBases(
            int refBaseExtensionDistance, int newMinAlignedPosition, int newMaxAlignedPosition, final RefBaseAssembly refBaseAssembly)
    {
        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int newBaseLength = mBases.length + refBaseExtensionDistance;

        boolean isForwardJunction = mJunction.isForward();

        if(isForwardJunction)
            mJunctionSequenceIndex += refBaseExtensionDistance;

        int baseOffset = isForwardJunction ? refBaseExtensionDistance : 0;

        int refBaseOffset = refBaseAssembly != null && isForwardJunction ? 0 : extensionLength();

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        mMinAlignedPosition = newMinAlignedPosition;
        mMaxAlignedPosition = newMaxAlignedPosition;

        // forward junction: 0 up to base offset - 1 -> set to zero, base offset to new length -> copy bases
        // reverse junction: 0 up to base offset - 1 -> copy bases, base offset to new length -> set to zero

        for(int i = 0; i < newBaseLength; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    if(refBaseAssembly != null && i < refBaseAssembly.bases().length)
                        mBases[i] = refBaseAssembly.bases()[i];
                    else
                        mBases[i] = 0;

                    mBaseQuals[i] = 0;
                }
                else
                {
                    mBases[i] = existingBases[i - baseOffset];
                    mBaseQuals[i] = existingQuals[i - baseOffset];
                }
            }
            else
            {
                if(i < existingBases.length)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                }
                else
                {
                    int refBaseIndex = i - refBaseOffset;
                    if(refBaseAssembly != null && refBaseIndex < refBaseAssembly.bases().length)
                        mBases[i] = refBaseAssembly.bases()[refBaseIndex];
                    else
                        mBases[i] = 0;

                    mBaseQuals[i] = 0;
                }
            }
        }
    }

    public void mergeRefBaseAssembly(final RefBaseAssembly refBaseAssembly)
    {
        // find the longest length of aligned reference bases extending back from the junction
        int newRefBaseCount = refBaseAssembly.validRefBaseLength();
        int existingRefBaseCount = refBaseLength();

        if(newRefBaseCount > existingRefBaseCount)
        {
            if(isForwardJunction())
                mMinAlignedPosition = refBaseAssembly.minAlignedPosition();
            else
                mMaxAlignedPosition = refBaseAssembly.maxAlignedPosition();

            int refBaseExtension = newRefBaseCount - existingRefBaseCount;

            extendBases(refBaseExtension, mMinAlignedPosition, mMaxAlignedPosition, refBaseAssembly);
        }

        for(SupportRead support : refBaseAssembly.support())
        {
            addRead(support.cachedRead(), support.type(), null);

            // once added, clear the cached data read
            support.clearCachedRead();
        }

        // use the ref assembly to fill in any missing bases

        List<int[]> emptyBaseRanges = findUnsetBases(mBases);

        if(!emptyBaseRanges.isEmpty())
        {
            SV_LOGGER.debug("assembly({}) refBases(existing={} new={}) empty ranges: {}",
                    refBaseAssembly, existingRefBaseCount, newRefBaseCount, emptyBaseRanges);
        }
    }

    private int[] getReadIndexCoordinates(final Read read, final int readJunctionIndex, boolean byReferenceBases)
    {
        // using the position of the junction within the read coordinates, gets the read index range either from the junction into
        // soft-clip / junction bases, or from the junction, but not including it, back into reference bases
        int[] readIndexCoords = {0, 0};

        if(byReferenceBases)
        {
            if(mJunction.isForward())
            {
                readIndexCoords[0] = 0;
                readIndexCoords[1] = readJunctionIndex - 1; // the base at the junction will have already been set
            }
            else
            {
                readIndexCoords[0] = readJunctionIndex + 1;
                readIndexCoords[1] = read.basesLength() - 1;
            }
        }
        else
        {
            if(mJunction.isForward())
            {
                readIndexCoords[0] = readJunctionIndex;
                readIndexCoords[1] = read.basesLength() - 1;
            }
            else
            {
                readIndexCoords[0] = 0;
                readIndexCoords[1] = readJunctionIndex;
            }
        }

        return readIndexCoords;
    }

    private int getReadAssemblyStartIndex(final int readJunctionIndex, final int readStartIndex, final boolean byReferenceBases)
    {
        if(byReferenceBases)
        {
            if(mJunction.isForward())
                return mJunctionSequenceIndex - (readJunctionIndex - readStartIndex);
            else
                return mJunctionSequenceIndex + 1;
        }
        else
        {
            if(mJunction.isForward())
                return mJunctionSequenceIndex; // ie bases starting from the junction
            else
                return mJunctionSequenceIndex - (readJunctionIndex - readStartIndex);
        }
    }

    public void removeSupportReads(final Set<String> readIds)
    {
        int index = 0;

        while(index < mSupport.size())
        {
            if(readIds.contains(mSupport.get(index).id()))
                mSupport.remove(index);
            else
                ++index;
        }
    }

    public boolean hasReadSupport(final Read read)
    {
        return read != null && mSupport.stream().anyMatch(x -> x.cachedRead() == read);
    }

    public void clearSupportCachedRead() { mSupport.forEach(x -> x.clearCachedRead()); }

    // caching repeat info needs careful consideration since any extension of ref bases invalidates the values,
    // at least for +ve orientation assemblies
    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
        {
            mRepeatInfo.addAll(repeats);
            mRefBasesRepeatedTrimmed = RepeatInfo.buildTrimmedRefBaseSequence(this, MIN_VARIANT_LENGTH);
        }
        else
        {
            mRefBasesRepeatedTrimmed = formRefBaseSequence(MIN_VARIANT_LENGTH);
        }

        mRefBaseTrimLength = calcTrimmedRefBaseLength(this, MIN_VARIANT_LENGTH);
    }

    public String refBasesRepeatedTrimmed() { return mRefBasesRepeatedTrimmed; }
    public int refBaseTrimLength() { return mRefBaseTrimLength; }

    public List<RefSideSoftClip> refSideSoftClips() { return mRefSideSoftClips; }

    public boolean checkAddRefSideSoftClip(final Read read)
    {
        return RefSideSoftClip.checkAddRefSideSoftClip(mRefSideSoftClips, mJunction, read);
    }

    public List<RemoteRegion> remoteRegions() { return mRemoteRegions; }
    public void addRemoteRegions(final List<RemoteRegion> regions) { mRemoteRegions.addAll(regions); }

    public PhaseGroup phaseGroup() { return mPhaseGroup; }
    public void setPhaseGroup(final PhaseGroup phaseGroup) { mPhaseGroup = phaseGroup; }

    public PhaseSet phaseSet()
    {
        // can only be a part of one and could add a reference but for now is retrieved
        return mPhaseGroup != null ? mPhaseGroup.findPhaseSet(this) : null;
    }

    public AssemblyOutcome outcome() { return mOutcome; }

    public void setOutcome(final AssemblyOutcome outcome)
    {
        if(mOutcome != REMOTE_REGION && mOutcome != LOCAL_REF) // persist classification for now
            mOutcome = outcome;
    }

    public AlignmentOutcome alignmentOutcome() { return mAlignmentOutcome; }
    public void setAlignmentOutcome(final AlignmentOutcome outcome) { mAlignmentOutcome = outcome; }

    public JunctionAssembly(
            final JunctionAssembly initialAssembly, final RefSideSoftClip refSideSoftClip,
            final List<SupportRead> initialSupport, final Set<String> excludedReadIds)
    {
        // build a junction assembly from an initial junction where the ref bases are truncated due to branching (likely short TI)
        mJunction = initialAssembly.junction();

        // copy the initial assembly's extension bases and ref bases up to the ref-side soft clip
        int extensionLength = initialAssembly.extensionLength();
        int newBaseLength = extensionLength + abs(mJunction.Position - refSideSoftClip.Position) + 1;
        int initialBaseLength = initialAssembly.baseLength();

        int baseLengthDiff = initialBaseLength - newBaseLength;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int assemblyIndexOffset = 0;

        if(mJunction.isForward())
        {
            assemblyIndexOffset = baseLengthDiff;
            mJunctionSequenceIndex = newBaseLength - extensionLength - 1;
            mMinAlignedPosition = refSideSoftClip.Position;
            mMaxAlignedPosition = initialAssembly.maxAlignedPosition();
        }
        else
        {
            mJunctionSequenceIndex = initialAssembly.junctionIndex();
            mMinAlignedPosition = initialAssembly.minAlignedPosition();
            mMaxAlignedPosition = refSideSoftClip.Position;
        }

        int assemblyIndex = assemblyIndexOffset;
        for(int i = 0; i < mBases.length; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
            {
                mBases[i] = 0;
                mBaseQuals[i] = 0;
            }
            else
            {
                if(assemblyIndex >= initialAssembly.bases().length)
                    break;

                mBases[i] = initialAssembly.bases()[assemblyIndex];
                mBaseQuals[i] = initialAssembly.baseQuals()[assemblyIndex];
            }
        }

        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();

        mRepeatInfo = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = initialAssembly.refBasesRepeatedTrimmed();
        mRefBaseTrimLength = initialAssembly.refBaseTrimLength();
        mRefSideSoftClips = Lists.newArrayList(refSideSoftClip);
        mRemoteRegions = Lists.newArrayList();
        mFilters = Sets.newHashSet();
        mMergedAssemblies = 0;
        mMismatchReadCount = 0;
        mPhaseGroup = null;

        mStats = new AssemblyStats();

        SupportRead initialRead = null;

        for(SupportRead support : initialSupport)
        {
            if(excludedReadIds.contains(support.id()))
                continue;

            mSupport.add(support);

            mStats.addRead(support, mJunction, null);

            if(support.id().equals(initialAssembly.initialReadId()))
                initialRead = support;
        }

        mInitialReadId = initialRead != null ? initialRead.id() : (!mSupport.isEmpty() ? mSupport.get(0).id() : "");
        mIndelCoords = initialAssembly.indelCoords();
        mOutcome = initialAssembly.outcome();
        mAlignmentOutcome = initialAssembly.alignmentOutcome();
    }

    public void addCandidateSupport(final Read read, final SupportType type)
    {
        mCandidateSupport.add(new SupportRead(read, type, 0, 0, 0));
        ++mStats.CandidateSupportCount;
    }

    public List<SupportRead> candidateSupport() { return mCandidateSupport; }
    public void clearCandidateSupport() { mCandidateSupport.clear(); }

    public Set<FilterType> filters() { return mFilters; }
    public boolean passing() { return mFilters.isEmpty(); }
    public void addFilter(final FilterType filterType) { mFilters.add(filterType); }

    public String toString()
    {
        return format("junc(%s) coords(extLen=%d refBasePos=%d len=%d) juncIndex(%d) support(%d) mismatches(%d)",
                mJunction, extensionLength(), refBasePosition(), baseLength(), mJunctionSequenceIndex,
                mSupport.size(), mMismatchReadCount);
    }

    public boolean hasUnsetBases() { return !findUnsetBases(mBases).isEmpty(); }

    public String formFullSequence() { return formJunctionSequence(refBaseLength()); }

    public String formJunctionSequence()
    {
        return formJunctionSequence(0);
    }

    public String formJunctionSequence(final int includeRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mJunction.isForward())
        {
            seqIndexStart = max(0, mJunctionSequenceIndex + 1 - includeRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionSequenceIndex - 1 + includeRefBaseCount, mBases.length - 1);
        }

        return formSequence(seqIndexStart, seqIndexEnd);
    }

    public String formRefBaseSequence() { return formRefBaseSequence(-1); }

    public String formRefBaseSequence(int maxRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mJunction.isForward())
        {
            seqIndexEnd = mJunctionSequenceIndex;

            if(maxRefBaseCount > 0)
                seqIndexStart = seqIndexEnd - maxRefBaseCount + 1;
            else
                seqIndexStart = 0;
        }
        else
        {
            seqIndexStart = mJunctionSequenceIndex;

            if(maxRefBaseCount > 0)
                seqIndexEnd = seqIndexStart + maxRefBaseCount - 1;
            else
                seqIndexEnd = mBases.length - 1;
        }

        return formSequence(seqIndexStart, seqIndexEnd);
    }

    public String formSequence(int seqIndexStart, int seqIndexEnd)
    {
        StringBuilder sb = new StringBuilder();

        for(int index = max(seqIndexStart, 0); index <= min(seqIndexEnd, mBases.length - 1); ++index)
        {
            sb.append((char)mBases[index]);
        }

        return sb.toString();
    }

    @VisibleForTesting
    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] quals, final int junctionIndex)
    {
        mJunction = junction;
        mInitialReadId = null;

        mJunctionSequenceIndex = junctionIndex;

        // pos = 20, index = 10, length = 21, min align = 10, max align = 30

        if(mJunction.isForward())
        {
            mMinAlignedPosition = mJunction.Position - junctionIndex;
            mMaxAlignedPosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }
        else
        {
            mMinAlignedPosition = mJunction.Position - junctionIndex;
            mMaxAlignedPosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }

        mBases = copyArray(bases);
        mBaseQuals = copyArray(quals);
        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = "";
        mRefBaseTrimLength = 0;
        mRemoteRegions = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mFilters = Sets.newHashSet();
        mMergedAssemblies = 0;
        mOutcome = UNSET;
        mAlignmentOutcome = NO_SET;
        mMismatchReadCount = 0;
        mStats = new AssemblyStats();
        mIndelCoords = null;
    }

    @VisibleForTesting
    public void addJunctionRead(final Read read, boolean registerMismatches)
    {
        SupportRead support = new SupportRead(read, JUNCTION, 0, 0, 0);
        mSupport.add(support);
        mStats.addRead(support, mJunction, read);
    }
}
