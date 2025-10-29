package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.calcTrimmedRefBaseLength;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.readQualFromJunction;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.convertedIndelCrossesJunction;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findInsertedBases;
import static com.hartwig.hmftools.esvee.common.CommonUtils.aboveMinQual;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.types.AssemblyOutcome.UNSET;
import static com.hartwig.hmftools.esvee.assembly.types.RepeatInfo.findRepeats;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.INDEL;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.RefReadParseState;
import com.hartwig.hmftools.esvee.assembly.RefBaseSeqBuilder;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.CigarElement;

public class JunctionAssembly
{
    private int mAssemblyId;
    private final Junction mJunction;

    private int mJunctionIndex; // position of the junction in the read bases

    // aligned position on the ref base side
    private int mRefBasePosition;
    private final List<CigarElement> mRefBaseCigarElements;

    private IndelCoords mIndelCoords;
    private boolean mHasLineSequence;

    private byte mBases[];
    private byte mBaseQuals[];

    private final List<SupportRead> mSupport;
    private final List<Read> mCandidateSupport;
    private final List<Read> mConcordantCandidates; // mates of junction reads past the junction
    private final List<Read> mUnmappedCandidates;
    private final List<RefSideSoftClip> mRefSideSoftClips;

    private final List<RepeatInfo> mRepeatInfo;
    private String mRefBasesRepeatedTrimmed;
    private int mRefBaseTrimLength;

    private final List<RemoteRegion> mRemoteRegions;

    private PhaseGroup mPhaseGroup;
    private AssemblyOutcome mOutcome;
    private String mAssemblyAlignmentInfo;

    // info only
    private final String mInitialReadId;
    private int mMergedAssemblies;
    private int mMismatchReadCount;
    private String mExtBaseBuildInfo;

    private final AssemblyStats mStats;

    public JunctionAssembly(
            final Junction junction, final byte[] bases, final byte[] baseQualities, final List<SupportRead> assemblySupport,
            final List<RepeatInfo> repeatInfo)
    {
        mJunction = junction;

        // set initial bounds from the assembly support
        Read maxJunctionBaseQualRead = null;
        int maxJunctionBaseQualTotal = 0;

        IndelCoords indelCoords = null;
        mRefBasePosition = junction.Position; // initialised to the same prior to extending ref bases
        mRefBaseCigarElements = Lists.newArrayList();

        for(SupportRead support : assemblySupport)
        {
            Read read = support.cachedRead();

            int junctionBaseQualTotal = readQualFromJunction(read, junction);

            if(junctionBaseQualTotal > maxJunctionBaseQualTotal)
            {
                maxJunctionBaseQualTotal = junctionBaseQualTotal;
                maxJunctionBaseQualRead = read;
            }

            if(read.indelCoords() != null && read.indelCoords().Length >= MIN_INDEL_LENGTH && indelCoords == null)
            {
                indelCoords = read.indelCoords();

                if(indelCoords.isInsert())
                    indelCoords.setInsertedBases(findInsertedBases(read));
            }
        }

        mIndelCoords = indelCoords;
        mHasLineSequence = false;

        mInitialReadId = maxJunctionBaseQualRead != null ? maxJunctionBaseQualRead.id() :
                (!assemblySupport.isEmpty() ? assemblySupport.get(0).id() : "null");

        mBases = bases;
        mBaseQuals = baseQualities;

        mJunctionIndex = junction.isForward() ? 0 : mBases.length - 1;

        mSupport = Lists.newArrayList(assemblySupport);
        mCandidateSupport = Lists.newArrayList();
        mConcordantCandidates = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();
        mRepeatInfo = repeatInfo;
        mRefSideSoftClips = Lists.newArrayList();
        mRemoteRegions = Lists.newArrayList();
        mMergedAssemblies = 0;
        mPhaseGroup = null;
        mOutcome = UNSET;
        mAssemblyAlignmentInfo = null;
        mMismatchReadCount = 0;
        mExtBaseBuildInfo = null;

        mStats = new AssemblyStats();
    }

    public void setId(int id) { mAssemblyId = id; }
    public int id() { return mAssemblyId; }

    public Junction junction() { return mJunction; }
    public boolean isForwardJunction() { return mJunction.isForward(); }
    public boolean isReverseJunction() { return mJunction.isReverse(); }

    public boolean indel() { return mJunction.indelBased(); }
    public boolean discordantOnly() { return mJunction.DiscordantOnly; }
    public IndelCoords indelCoords() { return mIndelCoords; }

    public int mergedAssemblyCount() { return mMergedAssemblies; }
    public void addMergedAssembly() { ++mMergedAssemblies; }

    public int junctionIndex() { return mJunctionIndex; };
    public void setJunctionIndex(int index) { mJunctionIndex = index; };

    // eg 21 bases, junction index at 10 (so 0-9 = 10 before, 11-20 = 10 after), note: doesn't count the junction base
    public int lowerDistanceFromJunction() { return mJunctionIndex; };
    public int upperDistanceFromJunction() { return mBases.length - mJunctionIndex - 1; };

    public int refBaseLength() { return (mJunction.isForward() ? lowerDistanceFromJunction() : upperDistanceFromJunction()) + 1; }

    // base count beyond the junction
    public int extensionLength() { return mJunction.isForward() ? upperDistanceFromJunction() : lowerDistanceFromJunction(); }

    public int refBasePosition() { return mRefBasePosition; }
    public String refBaseCigar() { return CigarUtils.cigarElementsToStr(mRefBaseCigarElements); }
    public int baseLength() { return mBases.length; }

    public byte[] bases() { return mBases; }
    public byte[] baseQuals() { return mBaseQuals; }

    public String initialReadId() { return mInitialReadId; }

    public boolean hasLineSequence() { return mHasLineSequence; }
    public void markLineSequence() { mHasLineSequence = true; }
    public void unmarkLineSequence() { mHasLineSequence = false; }

    public List<SupportRead> support() { return mSupport; }
    public int supportCount() { return mSupport.size(); }

    public void addCandidateSupport(final Read read)
    {
        mCandidateSupport.add(read);
        ++mStats.CandidateSupportCount;
    }

    public List<Read> candidateSupport() { return mCandidateSupport; }

    public void addConcordantCandidate(final Read read) { mConcordantCandidates.add(read); }
    public List<Read> concordantCandidates() { return mConcordantCandidates; }

    public void addUnmappedRead(final Read read)
    {
        mUnmappedCandidates.add(read);
        ++mStats.UnmappedReadCount;
    }

    public List<Read> unmappedReads() { return mUnmappedCandidates; }

    public void clearCandidateSupport()
    {
        mCandidateSupport.clear();
        mUnmappedCandidates.clear();
        mConcordantCandidates.clear();
    }

    public AssemblyStats stats() { return mStats; }

    public int mismatchReadCount() { return mMismatchReadCount; }
    public void addMismatchReadCount(int count) { mMismatchReadCount += count; }

    public void addRead(final Read read, final ReadAssemblyIndices readAssemblyIndices, final SupportType type)
    {
        if(readAssemblyIndices == ReadAssemblyIndices.INVALID_INDICES)
            return;

        int mismatchCount = 0;
        int highQualMatchCount = 0;
        int assemblyIndex = readAssemblyIndices.AssemblyIndexStart;

        for(int i = readAssemblyIndices.ReadIndexStart; i <= readAssemblyIndices.ReadIndexEnd; ++i, ++assemblyIndex)
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

                if(aboveMinQual(qual))
                    ++highQualMatchCount;
            }
            else
            {
                if(mBases[assemblyIndex] == base || belowMinQual(qual))
                {
                    if(qual > mBaseQuals[assemblyIndex])
                        mBaseQuals[assemblyIndex] = qual;

                    if(aboveMinQual(qual))
                        ++highQualMatchCount;
                }
                else if(belowMinQual(mBaseQuals[assemblyIndex]) || qual > mBaseQuals[assemblyIndex])
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

        int junctionReadStartDistance = readAssemblyIndices.junctionReadStartDistance(mJunctionIndex);
        addSupport(read, type, junctionReadStartDistance, highQualMatchCount, mismatchCount);
    }

    public void addSupport(
            final Read read, final SupportType type, int junctionReadStartDistance, int matches, int mismatches)
    {
        boolean isIndelCrossingJunction = convertedIndelCrossesJunction(mJunction, read);
        SupportType adjustedType = type == JUNCTION && isIndelCrossingJunction ? INDEL : type;
        SupportRead support = new SupportRead(read, adjustedType, junctionReadStartDistance, matches, mismatches);

        mSupport.add(support);
    }

    public void setRefBases(final RefBaseSeqBuilder refBaseSeqBuilder)
    {
        mRefBasePosition = refBaseSeqBuilder.refBasePosition();
        mRefBaseCigarElements.addAll(refBaseSeqBuilder.cigarElements());

        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        // build out the ref base sequence
        int refBaseExtension = refBaseSeqBuilder.refBaseLength() - 1; // since already includes the ref base at the junction
        int newBaseLength = mBases.length + refBaseExtension;
        boolean isForwardJunction = mJunction.isForward();

        int refBaseIndex = isForwardJunction ? 0 : 1;
        int baseOffset = isForwardJunction ? refBaseExtension : 0;

        if(isForwardJunction)
            mJunctionIndex += refBaseExtension;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        for(int i = 0; i < mBases.length; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    mBases[i] = refBaseSeqBuilder.bases()[refBaseIndex];
                    mBaseQuals[i] = refBaseSeqBuilder.baseQualities()[refBaseIndex];
                    ++refBaseIndex;
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
                    mBases[i] = refBaseSeqBuilder.bases()[refBaseIndex];
                    mBaseQuals[i] = refBaseSeqBuilder.baseQualities()[refBaseIndex];
                    ++refBaseIndex;
                }
            }
        }

        // update the support info for ref base mismatches - can rely on the support reads matching
        for(int i = 0; i < mSupport.size(); ++i)
        {
            SupportRead read = mSupport.get(i);
            RefReadParseState readState = refBaseSeqBuilder.reads().get(i);

            if(readState.isValid() && !readState.exceedsMaxMismatches())
            {
                read.setReferenceMismatches(readState.mismatches());
                checkAddRefSideSoftClip(read.cachedRead());
            }
        }
    }

    public void trimRefBasePosition(int newRefBasePosition)
    {
        if(isForwardJunction())
        {
            int trimLength = newRefBasePosition - mRefBasePosition;

            if(trimLength <= 0 || trimLength >= mBases.length - 1)
                return;

            mBases = subsetArray(mBases, trimLength, mBases.length - 1);
            mBaseQuals = subsetArray(mBaseQuals, trimLength, mBaseQuals.length - 1);
            mJunctionIndex -= trimLength;
        }
        else
        {
            int trimLength = mRefBasePosition - newRefBasePosition;

            if(trimLength <= 0 || trimLength >= mBases.length - 1)
                return;

            mBases = subsetArray(mBases, 0, mBases.length - 1 - trimLength);
            mBaseQuals = subsetArray(mBaseQuals, 0, mBaseQuals.length - 1 - trimLength);
        }

        mRefBasePosition = newRefBasePosition;

        // note that the ref base cigar is not adjusted since it is curently not extended from additional ref based reads either
    }

    public void extendRefBases(int newRefBasePosition, final RefGenomeInterface refGenome)
    {
        // extend the number of ref bases to accommodate new ref base information from existing or new reads
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int refBaseExtension = abs(newRefBasePosition - mRefBasePosition);

        int newBaseLength = mBases.length + refBaseExtension;

        int existingRefBasePosition = mRefBasePosition;

        mRefBasePosition = newRefBasePosition;

        boolean isForwardJunction = mJunction.isForward();

        if(isForwardJunction)
            mJunctionIndex += refBaseExtension;

        int baseOffset = isForwardJunction ? refBaseExtension : 0;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        byte[] refGenomeBases = null;

        if(refGenome != null)
        {
            int refBaseStart, refBaseEnd;

            if(mJunction.isForward())
            {
                refBaseStart = newRefBasePosition;
                refBaseEnd = existingRefBasePosition - 1;
            }
            else
            {
                refBaseStart = existingRefBasePosition + 1;
                refBaseEnd = newRefBasePosition;
            }

            refGenomeBases = refGenome.getBases(mJunction.Chromosome, refBaseStart, refBaseEnd);
        }

        // forward junction: 0 up to base offset - 1 -> set to zero or ref genome bases, base offset to new length -> copy bases
        // reverse junction: 0 up to base offset - 1 -> copy bases, base offset to new length -> set to zero or ref genome bases

        int refBaseIndex = 0;

        for(int i = 0; i < newBaseLength; ++i)
        {
            if(isForwardJunction)
            {
                if(i < baseOffset)
                {
                    if(refGenomeBases != null && refBaseIndex < refGenomeBases.length)
                        mBases[i] = refGenomeBases[refBaseIndex];
                    else
                        mBases[i] = 0;

                    ++refBaseIndex;
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
                    if(refGenomeBases != null && refBaseIndex < refGenomeBases.length)
                        mBases[i] = refGenomeBases[refBaseIndex];
                    else
                        mBases[i] = 0;

                    ++refBaseIndex;
                    mBaseQuals[i] = 0;
                }
            }
        }
    }

    public void expandExtensionBases(final byte[] extensionBases, final byte[] extensionBaseQuals, final List<SupportRead> supportReads)
    {
        // extend the same of ref bases to accommodate new ref base information from existing or new reads
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int newExtensionLength = extensionBases.length;
        int baseLengthChange = newExtensionLength - extensionLength();
        int existingBaseLength = mBases.length;
        int newBaseLength = existingBaseLength + baseLengthChange;

        if(mJunction.isReverse())
            mJunctionIndex += baseLengthChange;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int newExtBaseIndex = 0;
        for(int i = 0; i < newBaseLength; ++i)
        {
            if(mJunction.isForward())
            {
                if(i < existingBaseLength)
                {
                    mBases[i] = existingBases[i];
                    mBaseQuals[i] = existingQuals[i];
                }
                else
                {
                    mBases[i] = extensionBases[newExtBaseIndex];
                    mBaseQuals[i] = extensionBaseQuals[newExtBaseIndex];
                }

                if(i > mJunctionIndex) // start iterating through the new extension bases once past the junction index
                    ++newExtBaseIndex;
            }
            else
            {
                if(i < baseLengthChange)
                {
                    mBases[i] = extensionBases[newExtBaseIndex];
                    mBaseQuals[i] = extensionBaseQuals[newExtBaseIndex];
                    ++newExtBaseIndex;
                }
                else
                {
                    mBases[i] = existingBases[i - baseLengthChange];
                    mBaseQuals[i] = existingQuals[i - baseLengthChange];
                }
            }
        }

        for(SupportRead support : supportReads)
        {
            mSupport.add(support);
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

    public void trimRefBases(int newRefPosition)
    {
        // trims to the specified ref position, and assumes no indels are included
        byte[] existingBases = copyArray(mBases);
        byte[] existingQuals = copyArray(mBaseQuals);

        int existingRefPosition = refBasePosition();

        int newBaseLength = mBases.length - abs(existingRefPosition - newRefPosition);
        int baseLengthChange = mBases.length - newBaseLength;

        if(mJunction.isForward())
            mJunctionIndex -= baseLengthChange;

        mRefBasePosition = newRefPosition;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int baseOffset = mJunction.isForward() ? baseLengthChange : 0;

        for(int i = 0; i < newBaseLength; ++i)
        {
            mBases[i] = existingBases[i + baseOffset];
            mBaseQuals[i] = existingQuals[i + baseOffset];
        }

        buildRepeatInfo();
    }

    public boolean hasReadSupport(final Read read)
    {
        return read != null && mSupport.stream().anyMatch(x -> x.cachedRead() == read);
    }

    public void clearSupportCachedReads()
    {
        for(SupportRead read : mSupport)
        {
            if(read.cachedRead() != null)
            {
                // register stats info while still has access to raw read
                mStats.addRead(read, mJunction, read.cachedRead());
                read.clearCachedRead();
            }
        }
    }

    // caching repeat info needs careful consideration since any extension of ref bases invalidates the values,
    // at least for +ve orientation assemblies
    public List<RepeatInfo> repeatInfo() { return mRepeatInfo; }

    public void buildRepeatInfo()
    {
        mRepeatInfo.clear();
        List<RepeatInfo> repeats = findRepeats(mBases);
        if(!repeats.isEmpty())
        {
            mRepeatInfo.addAll(repeats);
            mRefBasesRepeatedTrimmed = RepeatInfo.buildTrimmedRefBaseSequence(this, MIN_VARIANT_LENGTH);
        }
        else
        {
            mRefBasesRepeatedTrimmed = formRefBaseSequence(MIN_VARIANT_LENGTH);
        }

        mRefBaseTrimLength = calcTrimmedRefBaseLength(this);
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

    public void setOutcome(final AssemblyOutcome outcome) { setOutcome(outcome, false); }

    public void setOutcome(final AssemblyOutcome outcome, boolean override)
    {
        if(!override && mOutcome.ordinal() <= outcome.ordinal()) // only override if a stronger type of link
            return;

        mOutcome = outcome;
    }

    public void setAssemblyAlignmentInfo(final String info) { mAssemblyAlignmentInfo = info; }
    public String assemblyAlignmentInfo() { return mAssemblyAlignmentInfo != null ? mAssemblyAlignmentInfo : mJunction.coords(); }

    public void setExtBaseBuildInfo(final String info) { mExtBaseBuildInfo = info; }
    public String extBaseBuildInfo() { return mExtBaseBuildInfo != null ? mExtBaseBuildInfo : ""; }

    public JunctionAssembly(
            final JunctionAssembly initialAssembly, final RefSideSoftClip refSideSoftClip, int refBaseLength,
            final List<SupportRead> initialSupport)
    {
        // build a junction assembly from an initial junction where the ref bases are truncated due to branching (likely short TI)
        mJunction = initialAssembly.junction();

        // copy the initial assembly's extension bases and ref bases up to the ref-side soft clip
        int extensionLength = initialAssembly.extensionLength();
        int newBaseLength = extensionLength + refBaseLength;
        int initialBaseLength = initialAssembly.baseLength();
        int baseLengthDiff = initialBaseLength - newBaseLength;

        mBases = new byte[newBaseLength];
        mBaseQuals = new byte[newBaseLength];

        int assemblyIndexOffset = 0;

        if(mJunction.isForward())
        {
            assemblyIndexOffset = baseLengthDiff;
            mJunctionIndex = newBaseLength - extensionLength - 1;
        }
        else
        {
            mJunctionIndex = initialAssembly.junctionIndex();
        }

        mRefBasePosition = refSideSoftClip.Position;
        mRefBaseCigarElements = initialAssembly.mRefBaseCigarElements;

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
        mConcordantCandidates = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();

        mRepeatInfo = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = "";
        mRefBaseTrimLength = 0;
        mRefSideSoftClips = Lists.newArrayList(refSideSoftClip);
        mRemoteRegions = Lists.newArrayList();
        mMergedAssemblies = 0;
        mMismatchReadCount = 0;
        mPhaseGroup = null;

        mStats = new AssemblyStats();

        SupportRead initialRead = null;

        for(SupportRead support : initialSupport)
        {
            mSupport.add(support);

            mStats.addRead(support, mJunction, support.cachedRead());

            if(support.id().equals(initialAssembly.initialReadId()))
                initialRead = support;
        }

        mInitialReadId = initialRead != null ? initialRead.id() : (!mSupport.isEmpty() ? mSupport.get(0).id() : "");
        mIndelCoords = initialAssembly.indelCoords();
        mOutcome = UNSET;
    }

    public String toString()
    {
        return format("junc(%s) coords(extLen=%d refLen=%d refBasePos=%d len=%d juncIndex=%d) support(%d) mismatches(%d)",
                mJunction.coordsTyped(), extensionLength(), refBaseLength(), refBasePosition(), baseLength(), mJunctionIndex,
                mSupport.size(), mMismatchReadCount);
    }

    public String formFullSequence() { return new String(mBases); }

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
            seqIndexStart = max(0, mJunctionIndex + 1 - includeRefBaseCount);
            seqIndexEnd = mBases.length - 1;
        }
        else
        {
            seqIndexStart = 0;
            seqIndexEnd = min(mJunctionIndex - 1 + includeRefBaseCount, mBases.length - 1);
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
            seqIndexEnd = mJunctionIndex;

            if(maxRefBaseCount > 0)
                seqIndexStart = seqIndexEnd - maxRefBaseCount + 1;
            else
                seqIndexStart = 0;
        }
        else
        {
            seqIndexStart = mJunctionIndex;

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

    // convenience
    public int minAlignedPosition() { return mJunction.isForward() ? refBasePosition() : mJunction.Position; }
    public int maxAlignedPosition() { return mJunction.isReverse() ? refBasePosition() : mJunction.Position; }

    @VisibleForTesting
    public JunctionAssembly(final Junction junction, final byte[] bases, final byte[] quals, final int junctionIndex)
    {
        mJunction = junction;
        mInitialReadId = null;

        mJunctionIndex = junctionIndex;

        // pos = 20, index = 10, length = 21, min align = 10, max align = 30

        if(mJunction.isForward())
        {
            mRefBasePosition = mJunction.Position - junctionIndex;
        }
        else
        {
            mRefBasePosition = mJunction.Position + (bases.length - junctionIndex) - 1;
        }

        mRefBaseCigarElements = Lists.newArrayList();
        mBases = copyArray(bases);
        mBaseQuals = copyArray(quals);
        mSupport = Lists.newArrayList();
        mCandidateSupport = Lists.newArrayList();
        mConcordantCandidates = Lists.newArrayList();
        mUnmappedCandidates = Lists.newArrayList();
        mRepeatInfo = Lists.newArrayList();
        mRefBasesRepeatedTrimmed = "";
        mRefBaseTrimLength = 0;
        mRemoteRegions = Lists.newArrayList();
        mRefSideSoftClips = Lists.newArrayList();
        mMergedAssemblies = 0;
        mOutcome = UNSET;
        mMismatchReadCount = 0;
        mStats = new AssemblyStats();
        mIndelCoords = null;
    }

    @VisibleForTesting
    public void addJunctionRead(final Read read)
    {
        // positive if on the lower side of the junction
        int junctionReadStartDistance = mJunction.Position - read.unclippedStart();

        // assume all extension bases match
        int extensionBases = mJunction.isForward() ? read.unclippedEnd() - mJunction.Position : mJunction.Position - read.unclippedStart();

        SupportRead support = new SupportRead(read, JUNCTION, junctionReadStartDistance, extensionBases, 0);
        support.setReferenceMismatches(0);
        mSupport.add(support);
    }

    @VisibleForTesting
    public void setIndelCoords(final IndelCoords coords) { mIndelCoords = coords; }
}
