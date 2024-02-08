package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.flipOrientation;
import static com.hartwig.hmftools.esvee.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.common.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.read.ReadUtils.copyArray;

import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

public class JunctionAssembly
{
    private int mAssemblyId;
    private final Junction mJunction;
    private final Read mInitialRead;

    private int mJunctionSequenceIndex; // position of the junction in the read bases
    private int mMinAlignedPosition;
    private int mMaxAlignedPosition;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<AssemblySupport> mSupport;
    private final List<RefSideSoftClip> mRefSideSoftClips;

    private final SequenceMismatches mSequenceMismatches;
    private final List<RepeatInfo> mRepeatInfo;

    // info only
    private int mMergedAssemblies;

    private final List<RemoteRegion> mRemoteRegions;
    private RefBaseAssembly mRefBaseAssembly;

    private PrimaryPhaseGroup mPrimaryPhaseGroup;

    public JunctionAssembly(
            final Junction initialJunction, final Read read, final int maxExtensionDistance,
            final int minAlignedPosition, final int maxAlignedPosition)
    {
        mJunction = initialJunction;
        mInitialRead = read;

        mJunctionSequenceIndex = initialJunction.isForward() ? 0 : maxExtensionDistance;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        int initialAssemblyLength = maxExtensionDistance + 1;
        mBases = new byte[initialAssemblyLength];

        mBaseQuals = new byte[initialAssemblyLength];
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mRefBaseAssembly = null;
        mMergedAssemblies = 0;
        mPrimaryPhaseGroup = null;

        addInitialRead(read);
    }

    public void setId(int id) { mAssemblyId = id; }
    public int id() { return mAssemblyId; }

    public Read initialRead() { return mInitialRead; }
    public Junction junction() { return mJunction; }
    public boolean isForwardJunction() { return mJunction.isForward(); }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionSequenceIndex; };

    // eg 21 bases, junction index at 10 (so 0-9 before, 11-20 after)
    public int lowerDistanceFromJunction() { return mJunctionSequenceIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionSequenceIndex - 1; };
    public int refBaseLength() { return mJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction(); }

    // base count beyond the junction
    public int extensionLength() { return mJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    // purely for informational purposes at this stage - since past the junction is just the soft-clip length, and on the
    // reference side doesn't take into account any INDELs
    public int minAlignedPosition() { return mMinAlignedPosition; }
    public int maxAlignedPosition() { return mMaxAlignedPosition; }
    public int baseLength() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public List<AssemblySupport> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public SequenceMismatches mismatches() { return mSequenceMismatches; }

    private void addInitialRead(final Read read)
    {
        for(int i = 0; i < mBases.length; ++i)
        {
            mBases[i] = 0;
            mBaseQuals[i] = 0;
        }

        addRead(read, false);
    }

    public boolean checkReadMatches(final Read read, int permittedMismatches)
    {
        int mismatchCount = 0;

        int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);
        int[] readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, false);
        int assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], false);

        if(assemblyIndex < 0)
        {
            // SV_LOGGER.debug("readCoords({}) invalid for assembly({}) read({})", readCoords, toString(), read.toString());
            return false;
        }

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length) // CHECK: similar to issue with INDELs in addRead() ??
                break;

            if(!basesMatch(
                    read.getBases()[i], mBases[assemblyIndex], read.getBaseQuality()[i], mBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                ++mismatchCount;

                if(mismatchCount > permittedMismatches)
                    return false;
            }
        }

        return true;
   }

    public void addRead(final Read read, boolean registerMismatches)
    {
        addRead(read, registerMismatches, null);
    }

    public void extendReadSupport(final Read read, final AssemblySupport existingSupport)
    {
        addRead(read, true, existingSupport);
    }

    private void addRead(final Read read, boolean registerMismatches, @Nullable final AssemblySupport existingSupport)
    {
        int mismatchCount = 0;

        int readJunctionIndex;

        boolean byReference = existingSupport != null;

        if(existingSupport == null)
        {
            readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);
        }
        else
        {
            readJunctionIndex = existingSupport.junctionReadIndex();
        }

        int[] readIndexRange = getReadIndexCoordinates(read, readJunctionIndex, byReference);

        int assemblyIndex = getReadAssemblyStartIndex(readJunctionIndex, readIndexRange[SE_START], byReference);

        for(int i = readIndexRange[SE_START]; i <= readIndexRange[SE_END]; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= mBases.length)
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
            }
            else
            {
                if(mBases[assemblyIndex] == base || qual < LOW_BASE_QUAL_THRESHOLD)
                {
                    if((int)qual > (int)mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;
                }
                else if(mBaseQuals[assemblyIndex] < LOW_BASE_QUAL_THRESHOLD)
                {
                    mBases[assemblyIndex] = base;
                    mBaseQuals[assemblyIndex] = qual;
                }
                else
                {
                    ++mismatchCount;

                    if(registerMismatches)
                        mSequenceMismatches.add(assemblyIndex, base, read, qual);
                }
            }
        }

        if(existingSupport == null)
        {
            mSupport.add(new AssemblySupport(read, JUNCTION, assemblyIndex, readJunctionIndex, readIndexRange, mismatchCount));
        }
        else
        {
            int[] existingReadRange = existingSupport.readIndexRange();

            existingSupport.setReadIndexRange(
                    min(existingReadRange[SE_START], readIndexRange[SE_START]),
                    max(existingReadRange[SE_END], readIndexRange[SE_END]));

            existingSupport.setReferenceMismatches(mismatchCount);
        }
    }

    public void extendBases(int maxDistanceFromJunction, int newMinAlignedPosition, int newMaxAlignedPosition)
    {
        // copy existing bases and quals
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int newBaseLength = mBases.length + maxDistanceFromJunction;

        boolean isForwardJunction = mJunction.isForward();

        if(isForwardJunction)
            mJunctionSequenceIndex += maxDistanceFromJunction;

        int baseOffset = isForwardJunction ? maxDistanceFromJunction : 0;

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
                    mBases[i] = 0;
                    mBaseQuals[i] = 0;
                }
            }
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

    public void removeSupportRead(final Read read)
    {
        for(int i = 0; i < mSupport.size(); ++i)
        {
            if(mSupport.get(i).read() == read)
            {
                mSupport.remove(i);
                return;
            }
        }
    }

    public void checkAddReadSupport(final JunctionAssembly other)
    {
        for(AssemblySupport support : other.support())
        {
            if(hasReadSupport(support.read()))
                continue;

            addRead(support.read(), false);
        }
    }

    public boolean hasReadSupport(final Read read)
    {
        return read != null && mSupport.stream().anyMatch(x -> x.read() == read);
    }

    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(repeats != null)
            mRepeatInfo.addAll(repeats);
    }

    public List<RefSideSoftClip> refSideSoftClips() { return mRefSideSoftClips; }

    public boolean checkAddRefSideSoftClip(final Read read)
    {
        int refSideSoftClipPosition = -1;
        int refSideSoftClipLength = 0;

        if(isForwardJunction())
        {
            refSideSoftClipPosition = read.alignmentStart();
            refSideSoftClipLength = read.leftClipLength();
        }
        else
        {
            refSideSoftClipPosition = read.alignmentEnd();
            refSideSoftClipLength = read.rightClipLength();
        }

        if(refSideSoftClipLength == 0)
            return false;

        int softClipPosition = refSideSoftClipPosition;
        RefSideSoftClip existing = mRefSideSoftClips.stream().filter(x -> x.Position == softClipPosition).findFirst().orElse(null);

        if(existing == null)
        {
            mRefSideSoftClips.add(new RefSideSoftClip(softClipPosition, flipOrientation(mJunction.Orientation), read, refSideSoftClipLength));
        }
        else
        {
            existing.addRead(read, refSideSoftClipLength);
        }

        return true;
    }

    public void purgeRefSideSoftClips(int minCount, int minLength, int nonSoftClipRefPosition)
    {
        if(mRefSideSoftClips.isEmpty())
            return;

        // drop if not sufficient support or matches the original assembly's ref extension position anyway
        // or is close to it - where say homology causes aligned bases to match the soft-clip
        int index = 0;
        RefSideSoftClip matching = null;

        while(index < mRefSideSoftClips.size())
        {
            RefSideSoftClip refSideSoftClip = mRefSideSoftClips.get(index);

            boolean isPositionMatchOrClose = abs(refSideSoftClip.Position - nonSoftClipRefPosition) < 4;

            if(refSideSoftClip.readCount() < minCount || refSideSoftClip.maxLength() < minLength || isPositionMatchOrClose)
            {
                mRefSideSoftClips.remove(index);

                if(isPositionMatchOrClose)
                {
                    if(refSideSoftClip.Position == nonSoftClipRefPosition)
                        matching = refSideSoftClip;
                    else if(matching == null || matching.readCount() < refSideSoftClip.readCount())
                         matching = refSideSoftClip;
                }
            }
            else
            {
                ++index;
            }
        }

        // retain info about any matching soft-clip
        if(matching != null)
        {
            matching.markMatchesOriginal();
            mRefSideSoftClips.add(matching);
        }
    }

    public List<RemoteRegion> remoteRegions() { return mRemoteRegions; }
    public void addRemoteRegions(final List<RemoteRegion> regions) { mRemoteRegions.addAll(regions); }

    public RefBaseAssembly refBaseAssembly() { return mRefBaseAssembly; }
    public void setRefBaseAssembly(final RefBaseAssembly refBaseAssembly) { mRefBaseAssembly = refBaseAssembly; }

    public PrimaryPhaseGroup primaryPhaseGroup() { return mPrimaryPhaseGroup; }
    public void setPrimaryPhaseGroup(final PrimaryPhaseGroup primaryPhaseGroup) { mPrimaryPhaseGroup = primaryPhaseGroup; }

    public JunctionAssembly(
            final JunctionAssembly initialAssembly, final RefSideSoftClip refSideSoftClip, final Set<Read> excludedReads)
    {
        mJunction = initialAssembly.junction();
        mInitialRead = initialAssembly.initialRead(); // possible that this read isn't in the final assembly if it supports a different ref

        mJunctionSequenceIndex = initialAssembly.junctionIndex();

        mMinAlignedPosition = initialAssembly.minAlignedPosition();
        mMaxAlignedPosition = initialAssembly.maxAlignedPosition();

        mBases = new byte[initialAssembly.baseLength()];
        mBaseQuals = new byte[initialAssembly.baseLength()];

        mBases = copyArray(initialAssembly.bases());
        mBaseQuals = copyArray(initialAssembly.baseQuals());
        mSequenceMismatches = new SequenceMismatches();

        mSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList(refSideSoftClip);
        mRemoteRegions = Lists.newArrayList();
        mRefBaseAssembly = null;
        mMergedAssemblies = 0;
        mPrimaryPhaseGroup = null;

        for(AssemblySupport support : initialAssembly.support())
        {
            if(!excludedReads.contains(support.read()))
                mSupport.add(support);
        }
    }

    public String toString()
    {
        return format("junc(%s) range(%d - %d len=%d) juncIndex(%d) support(%d) mismatches(pos=%d all=%d)",
                mJunction, mMinAlignedPosition, mMaxAlignedPosition, baseLength(), mJunctionSequenceIndex,
                mSupport.size(), mSequenceMismatches.positionCount(), mSequenceMismatches.distinctBaseCount());
    }

    public String formSequence(final int maxRefBaseCount)
    {
        int seqIndexStart;
        int seqIndexEnd;

        if(mJunction.isForward())
        {
            seqIndexStart = max(0, mJunctionSequenceIndex - maxRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionSequenceIndex + maxRefBaseCount, mBases.length - 1);
        }

        StringBuilder sb = new StringBuilder();

        for(int index = seqIndexStart; index <= seqIndexEnd; ++index)
        {
            sb.append((char)mBases[index]);
        }

        return sb.toString();
    }


    @VisibleForTesting
    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] quals, final int minAlignedPosition, final int maxAlignedPosition)
    {
        mJunction = junction;
        mInitialRead = null;

        mJunctionSequenceIndex = junction.isForward() ? junction.Position - minAlignedPosition : maxAlignedPosition - junction.Position;

        mMinAlignedPosition = minAlignedPosition;
        mMaxAlignedPosition = maxAlignedPosition;

        mBases = copyArray(bases);
        mBaseQuals = copyArray(quals);
        mSequenceMismatches = new SequenceMismatches();
        mSupport = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mMergedAssemblies = 0;
    }
}
