package com.hartwig.hmftools.redux.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.bam.UmiReadType.DUAL;
import static com.hartwig.hmftools.common.bam.UmiReadType.SINGLE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentUtils.readToString;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroup
{
    private final String mUmiId; // the UMI if enabled
    private final List<Fragment> mFragments;
    private List<String> mReadIds;
    private int mFragmentCount;

    // reads from each fragment are organised into their like-types from which consensus reads can be formed
    private final List<SAMRecord>[] mReadGroups;
    private final boolean[] mReadGroupComplete;
    private final ReadTypeId[] mPrimaryReadTypeIndex; // details for primary and mate reads
    private final String mCoordinatesKey;

    private SAMRecord mPrimaryTemplateRead; // read on which the primary consensus read is based
    private String mGroupReadId;
    private boolean mDualStrand;

    private static final int MAX_READ_TYPES = ReadType.values().length;
    private static final int PRIMARY_READ_TYPES = ReadType.MATE.ordinal() + 1;

    public DuplicateGroup(final String id, final Fragment fragment)
    {
        mUmiId = id;
        mFragments = Lists.newArrayList(fragment);
        mReadIds = null;
        mReadGroups = new List[MAX_READ_TYPES];
        mReadGroupComplete = new boolean[MAX_READ_TYPES];
        mPrimaryReadTypeIndex = new ReadTypeId[PRIMARY_READ_TYPES];
        mFragmentCount = 0;
        mCoordinatesKey = fragment.coordinates().keyOriented();
        mPrimaryTemplateRead = null;
        mGroupReadId = null;
        mDualStrand = false;
    }

    public List<Fragment> fragments() { return mFragments; }
    public void addFragment(final Fragment fragment) { mFragments.add(fragment); }
    public int fragmentCount() { return mFragmentCount > 0 ? mFragmentCount : mFragments.size(); }

    public String coordinatesKey() { return mCoordinatesKey; }
    public FragmentCoordinates fragmentCoordinates() { return !mFragments.isEmpty() ? mFragments.get(0).coordinates() : null; }

    public String umiId() { return mUmiId; }

    public void registerDualStrand() { mDualStrand = true; }
    public boolean hasDualStrand() { return mDualStrand; }

    public void categoriseReads()
    {
        if(mFragmentCount > 0)
            return;

        mFragmentCount = mFragments.size();
        mReadIds = Lists.newArrayListWithExpectedSize(mFragmentCount);

        // establish expected lists by read type
        Fragment firstFragment = mFragments.get(0);

        mReadGroups[ReadType.PRIMARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

        for(int i = 0; i < firstFragment.reads().size(); ++i)
        {
            SAMRecord read = firstFragment.reads().get(i);

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            boolean hasValidSupp = suppData != null && HumanChromosome.contains(suppData.Chromosome);

            if(i == 0)
            {
                if(read.getReadPairedFlag() && HumanChromosome.contains(read.getMateReferenceName()))
                    mReadGroups[ReadType.MATE.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);

                if(hasValidSupp)
                    mReadGroups[ReadType.PRIMARY_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }
            else if(!read.getSupplementaryAlignmentFlag() && hasValidSupp)
            {
                mReadGroups[ReadType.MATE_SUPPLEMENTARY.ordinal()] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }
        }

        for(Fragment fragment : mFragments)
        {
            mReadIds.add(fragment.id());

            // FIXME: at the moment this fragment ID used to both store the UMI, and to indicate at the fragment is part of a duplicate group,
            // including for when UMIs are disabled and it's a standard duplicate group. Consider renaming or altering meaning
            fragment.setUmi(mUmiId != null ? mUmiId : "");
            fragment.setStatus(DUPLICATE);

            boolean checkFirstInPair = fragment.isPreciseInversion();

            // add non-supps first to establish the correct primary read type info
            fragment.reads().stream().filter(x -> !x.getSupplementaryAlignmentFlag()).forEach(x -> addRead(x, checkFirstInPair));
        }

        // add supplementaries once all primaries have been added and their expected supps & types registered
        mFragments.forEach(x -> x.reads().stream().filter(y -> y.getSupplementaryAlignmentFlag()).forEach(y -> addRead(y)));

        mFragments.clear();
    }

    public void addRead(final SAMRecord read) { addRead(read, false); }

    private void addRead(final SAMRecord read, boolean checkFirstInPair)
    {
        int readTypeIndex = getReadTypeIndex(read, checkFirstInPair);

        if(mReadGroups[readTypeIndex] == null)
        {
            mReadGroups[readTypeIndex] = Lists.newArrayListWithExpectedSize(mFragmentCount);
        }

        mReadGroups[readTypeIndex].add(read);
    }

    private enum ReadType
    {
        PRIMARY,
        MATE,
        PRIMARY_SUPPLEMENTARY,
        MATE_SUPPLEMENTARY;
    }

    private class ReadTypeId
    {
        // a class for distinguishing between reads in a fragment - first & second, and any supplementaries
        // first/second in pair is not sufficient where reverse orientation fragments have R1 and R2 swapped
        // so coordinates are also used
        public final String Chromosome;
        public final int UnclippedPosition;
        public final int[] PositionRange;
        public final byte Orientation;
        public final boolean FirstInPair;
        public boolean HasSupplementary;
        public final boolean Unmapped;

        public ReadTypeId(
                final String chromosome, final int unclippedPosition, final int position, final byte orientation,
                final boolean hasSupplementary, final boolean firstInPair, final boolean unmapped)
        {
            Chromosome = chromosome;
            UnclippedPosition = unclippedPosition;
            PositionRange = new int[] { position, position };
            Orientation = orientation;
            HasSupplementary = hasSupplementary;
            FirstInPair = firstInPair;
            Unmapped = unmapped;
        }

        public void updatePosition(final int position)
        {
            PositionRange[SE_START] = min(PositionRange[SE_START], position);
            PositionRange[SE_END] = max(PositionRange[SE_END], position);
        }

        private boolean primaryMatches(final SAMRecord read, boolean checkFirstInPair)
        {
            if(Unmapped || read.getReadUnmappedFlag())
                return Unmapped == read.getReadUnmappedFlag();

            if(!read.getReferenceName().equals(Chromosome))
                return false;

            if(orientation(read) != Orientation)
                return false;

            if(checkFirstInPair && FirstInPair != read.getFirstOfPairFlag())
                return false;

            return getFivePrimeUnclippedPosition(read) == UnclippedPosition;
        }

        public boolean supplementaryMatches(final SAMRecord read, final SupplementaryReadData suppData)
        {
            if(!suppData.Chromosome.equals(Chromosome))
                return false;

            if(orientation(read) != Orientation)
                return false;

            return positionWithin(suppData.Position, PositionRange[SE_START], PositionRange[SE_END]);
        }

        public String toString()
        {
            return format("%s:%d%s%s %d-%d %s",
                Chromosome, UnclippedPosition, Orientation == POS_ORIENT ? " " : "_R ",
                    FirstInPair ? "R1" : "R2", PositionRange[SE_START], PositionRange[SE_END],
                    HasSupplementary ? "has-supp" : "");
        }
    }

    private int getReadTypeIndex(final SAMRecord read, boolean checkFirstInPair)
    {
        if(!read.getSupplementaryAlignmentFlag())
        {
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);
            boolean hasValidSupp = suppData != null && HumanChromosome.contains(suppData.Chromosome);

            int index = 0;
            for(; index < mPrimaryReadTypeIndex.length; ++index)
            {
                if(mPrimaryReadTypeIndex[index] == null)
                    break;

                if(mPrimaryReadTypeIndex[index].primaryMatches(read, checkFirstInPair))
                {
                    // update details for this primary index
                    mPrimaryReadTypeIndex[index].updatePosition(read.getAlignmentStart());
                    mPrimaryReadTypeIndex[index].HasSupplementary |= hasValidSupp;
                    return index;
                }
            }

            if(index >= mPrimaryReadTypeIndex.length)
            {
                RD_LOGGER.error("group({}) non-supp read unmatched: {}", this, readToString(read));
                return index;
            }

            mPrimaryReadTypeIndex[index] = new ReadTypeId(
                    read.getReferenceName(),
                    read.getReadUnmappedFlag() ? 0 : getFivePrimeUnclippedPosition(read), read.getAlignmentStart(),
                    read.getReadUnmappedFlag() ? 0 : orientation(read),
                    hasValidSupp, read.getFirstOfPairFlag(), read.getReadUnmappedFlag());

            // register the expected supplementary read group
            if(hasValidSupp)
            {
                int suppReadTypeIndex = index == ReadType.PRIMARY.ordinal() ?
                        ReadType.PRIMARY_SUPPLEMENTARY.ordinal() : ReadType.MATE_SUPPLEMENTARY.ordinal();

                if(mReadGroups[suppReadTypeIndex] == null)
                    mReadGroups[suppReadTypeIndex] = Lists.newArrayListWithExpectedSize(mFragmentCount);
            }

            return index;
        }

        SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

        int matchedPrimaryIndex;

        if(suppData == null)
        {
            // logical assert since all supp should have this attribute set
            matchedPrimaryIndex = mPrimaryReadTypeIndex[0] != null && mPrimaryReadTypeIndex[0].HasSupplementary ? 0 : 1;
        }
        else
        {
            // match the supplementary details to that of the primary read
            matchedPrimaryIndex = mPrimaryReadTypeIndex[0] != null && mPrimaryReadTypeIndex[0].supplementaryMatches(read, suppData) ? 0 : 1;
        }

        return matchedPrimaryIndex == 0 ? ReadType.PRIMARY_SUPPLEMENTARY.ordinal() : ReadType.MATE_SUPPLEMENTARY.ordinal();
    }

    public boolean allReadsReceived()
    {
        for(int i = 0; i < mReadGroups.length; ++i)
        {
            if(mReadGroupComplete[i])
                continue;

            List<SAMRecord> readGroup = mReadGroups[i];

            if(readGroup != null && readGroup.size() < mFragmentCount)
                return false;
        }

        return true;
    }

    public boolean hasCompleteReadGroup()
    {
        return Arrays.stream(mReadGroups).anyMatch(x -> x != null && !x.isEmpty() && x.size() >= mFragmentCount);
    }

    public List<SAMRecord> popCompletedReads(final ConsensusReads consensusReads, boolean processIncompletes)
    {
        // take each read group type in turn and if complete, or in a final processing step, create a consensus read
        // and collect all reads to be written and then cleared from this UMI group
        List<SAMRecord> reads = Lists.newArrayList();

        for(int i = 0; i < mReadGroups.length; ++i)
        {
            if(mReadGroupComplete[i] && !processIncompletes)
                continue;

            List<SAMRecord> readGroup = mReadGroups[i];

            if(readGroup == null || readGroup.isEmpty())
                continue;

            if(readGroup.size() < mFragmentCount && !processIncompletes)
                continue;

            reads.addAll(readGroup);

            try
            {
                ConsensusReadInfo consensusReadInfo;

                if(i == ReadType.PRIMARY_SUPPLEMENTARY.ordinal() || i == ReadType.MATE_SUPPLEMENTARY.ordinal())
                {
                    // supplementaries can go to difference places and some reads have more than one, so go with the most frequent
                    consensusReadInfo = consensusReads.createConsensusRead(
                            findConsistentSupplementaries(readGroup), mPrimaryTemplateRead, mGroupReadId, mUmiId);
                }
                else
                {
                    consensusReadInfo = consensusReads.createConsensusRead(readGroup, mPrimaryTemplateRead, mGroupReadId, mUmiId);

                    // cache to ensure subsequent consensus reads use the same properties
                    if(mPrimaryTemplateRead == null)
                    {
                        mPrimaryTemplateRead = consensusReadInfo.TemplateRead;
                        mGroupReadId = consensusReadInfo.ConsensusRead.getReadName();

                    }
                }

                // set consensus read attributes
                int firstInPairCount = (int)readGroup.stream().filter(x -> x.getFirstOfPairFlag()).count();
                int readCount = readGroup.size();
                boolean isPrimaryGroup = (i == ReadType.PRIMARY.ordinal() || i == ReadType.PRIMARY_SUPPLEMENTARY.ordinal());

                if(!isPrimaryGroup)
                    firstInPairCount = readCount - firstInPairCount; // adjusted so both reads report the same ratio

                UmiReadType umiReadType = mDualStrand ? DUAL : SINGLE;

                addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount, firstInPairCount,  umiReadType);

                reads.add(consensusReadInfo.ConsensusRead);
            }
            catch(Exception e)
            {
                RD_LOGGER.error("error forming consensus: {}", toString());

                for(SAMRecord read : readGroup)
                {
                    RD_LOGGER.error("read: {}", readToString(read));
                }

                e.printStackTrace();
            }

            readGroup.clear();
            mReadGroupComplete[i] = true;
        }

        return reads;
    }

    private List<SAMRecord> findConsistentSupplementaries(final List<SAMRecord> readGroup)
    {
        Map<String,Integer> posCounts = Maps.newHashMap();

        String[] readChrPositions = new String[readGroup.size()];

        for(int i = 0; i < readGroup.size(); ++i)
        {
            SAMRecord read = readGroup.get(i);
            String chrPosStr = format("%s_%d", read.getReferenceName(), getFivePrimeUnclippedPosition(read));
            readChrPositions[i] = chrPosStr;
            Integer count = posCounts.get(chrPosStr);
            posCounts.put(chrPosStr, count != null ? count + 1 : 1);
        }

        if(posCounts.size() == 1)
            return readGroup;

        if(posCounts.size() == 2)
            return Lists.newArrayList(readGroup.get(0));

        String maxChrPos = "";
        int maxCount = 0;
        for(Map.Entry<String,Integer> entry : posCounts.entrySet())
        {
            if(entry.getValue() > maxCount)
            {
                maxChrPos = entry.getKey();
                maxCount = entry.getValue();
            }
        }

        List<SAMRecord> selectedReads = Lists.newArrayListWithExpectedSize(maxCount);
        for(int i = 0; i < readGroup.size(); ++i)
        {
            if(readChrPositions[i].equals(maxChrPos))
                selectedReads.add(readGroup.get(i));
        }

        return selectedReads;
    }

    public List<String> getReadIds()
    {
        if(!mFragments.isEmpty())
            return mFragments.stream().map(x -> x.id()).collect(Collectors.toList());
        else
            return mReadIds;
    }

    public int cachedReadCount()
    {
        return Arrays.stream(mReadGroups).filter(x -> x != null).mapToInt(x -> x.size()).sum();
    }

    public String toString()
    {
        if(mFragmentCount == 0)
            return format("id(%s) fragments(%d) coords(%s)", mUmiId, mFragments.size(), mCoordinatesKey);

        StringJoiner sj = new StringJoiner(", ");
        for(ReadType readType : ReadType.values())
        {
            int readTypeIndex = readType.ordinal();
            List<SAMRecord> readGroup = mReadGroups[readTypeIndex];

            if(readGroup == null)
                continue;

            String state = mReadGroupComplete[readTypeIndex] ?
                    (readGroup.isEmpty() ? "complete" : "extras") : (readGroup.size() == mFragmentCount ? "ready" : "pending");

            sj.add(format("%s=%d %s", readType, readGroup.size(), state));
        }

        return format("id(%s) fragments(%d) coords(%s) readCounts(%s)", mUmiId, mFragmentCount, mCoordinatesKey, sj);
    }
}
