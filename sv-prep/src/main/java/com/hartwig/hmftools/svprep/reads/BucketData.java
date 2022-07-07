package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.subtractExact;

import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.maxDeleteLength;
import static com.hartwig.hmftools.svprep.SvConstants.JUNCTION_SUPPORT_CAP;
import static com.hartwig.hmftools.svprep.SvConstants.LOW_BASE_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.SUPPORTING_READ_DISTANCE;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SoftClipSide;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.HotspotCache;

import htsjdk.samtools.CigarElement;

public class BucketData
{
    private final ChrBaseRegion mRegion;
    private final int mId;

    private final Map<String,ReadGroup> mReadGroups; // keyed by readId to additional reads (often supps) can be added
    private final List<ReadRecord> mSupportingReads; // for now only store as read

    private final List<JunctionData> mJunctions;
    private int mInitialSupportingReadCount;

    public BucketData(final int id, final ChrBaseRegion region)
    {
        mRegion = region;
        mId = id;
        mReadGroups = Maps.newHashMap();
        mSupportingReads = Lists.newArrayList();
        mJunctions = Lists.newArrayList();
        mInitialSupportingReadCount = 0;
    }

    public int id() { return mId; }
    public ChrBaseRegion region() { return mRegion; }

    public Collection<ReadGroup> readGroups() { return mReadGroups.values(); }
    public List<ReadRecord> supportingReads() { return mSupportingReads; }
    public List<JunctionData> junctions() { return mJunctions; }
    public int initialSupportingReadCount() { return mInitialSupportingReadCount; }

    public void addReadGroup(final ReadGroup readGroup)
    {
        ReadGroup existingGroup = mReadGroups.get(readGroup.id());

        if(existingGroup != null)
        {
            readGroup.reads().forEach(x -> existingGroup.addRead(x));
            return;
        }

        mReadGroups.put(readGroup.id(), readGroup);
    }

    public void addSupportingRead(final ReadRecord read)
    {
        mSupportingReads.add(read); // may be reassigned to a read group when the bucket is processed
    }

    public void filterJunctions(final HotspotCache hotspotCache, final ReadFilterConfig filterConfig)
    {
        int index = 0;
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.totalSupport() >= filterConfig.MinJunctionSupport)
            {
                ++index;
                continue;
            }

            // check for a hotspot match
            boolean matchesHotspot = junctionData.RemoteJunctions.stream()
                    .anyMatch(x -> hotspotCache.matchesHotspot(
                            mRegion.Chromosome, x.Chromosome, junctionData.Position, x.Position, junctionData.Orientation, x.Orientation));

            if(matchesHotspot && junctionData.totalSupport() >= MIN_HOTSPOT_JUNCTION_SUPPORT)
            {
                junctionData.markHotspot();
                ++index;
                continue;
            }

            mJunctions.remove(index);
        }
    }

    public void assignJunctionReads(final ReadFilterConfig filterConfig)
    {
        // first add any filtered reads to passing read groups
        // attempt to assign reads which support a junction, otherwise discard to keep to move to next bucket
        int index = 0;
        while(index < mSupportingReads.size())
        {
            ReadRecord read = mSupportingReads.get(index);
            ReadGroup existingGroup = mReadGroups.get(read.id());

            if(existingGroup != null)
            {
                existingGroup.addRead(read);
                mSupportingReads.remove(index);
            }
            else
            {
                ++index;
            }
        }

        // create junctions from the reads that passed the original filters
        createJunctions(filterConfig);

        // then find reads that support those junctions
        assignSupportingReads(filterConfig);
    }

    private void createJunctions(final ReadFilterConfig filterConfig)
    {
        // convert
        for(ReadGroup readGroup : mReadGroups.values())
        {
            // ignore any group with a short overlapping fragment, likely adapter
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), INSERT_MAP_OVERLAP)))
                continue;

            List<JunctionData> junctions = Lists.newArrayList();
            RemoteJunction remoteJunction = null;
            ReadRecord remoteJunctionRead = null;

            for(ReadRecord read : readGroup.reads())
            {
                if(read.hasSuppAlignment())
                    remoteJunction = RemoteJunction.fromSupplementaryData(read.supplementaryAlignment());

                handleInternalDelete(readGroup, read, filterConfig.MinDeleteLength);

                SoftClipSide scSide = SoftClipSide.fromCigar(read.cigar());

                if(scSide == null || scSide.Length < filterConfig.MinSoftClipLength)
                    continue;

                byte orientation = scSide.isLeft() ? NEG_ORIENT : POS_ORIENT;
                int position = scSide.isLeft() ? read.start() : read.end();

                if(!mRegion.containsPosition(position))
                {
                    // will only be cached if no junction is within this region
                    remoteJunction = new RemoteJunction(mRegion.Chromosome, position, orientation);
                    remoteJunctionRead = read;
                }
                else
                {
                    JunctionData junctionData = getOrCreateJunction(read, orientation);

                    if(!junctions.contains(junctionData))
                        junctions.add(junctionData);
                }
            }

            if(junctions.isEmpty() && remoteJunction != null && remoteJunction.Chromosome.equals(mRegion.Chromosome))
            {
                // convert this junction in a latter bucket to a junction for this group
                junctions.add(getOrCreateJunction(remoteJunctionRead, remoteJunction.Orientation));
                remoteJunction = null;
            }

            if(junctions.isEmpty())
                continue;

            junctions.forEach(x -> x.JunctionGroups.add(readGroup));

            if(remoteJunction != null)
            {
                final RemoteJunction newRemote = remoteJunction;

                for(JunctionData junctionData : junctions)
                {
                    if(junctionData.RemoteJunctions.stream().noneMatch(x -> x.matches(newRemote)))
                        junctionData.RemoteJunctions.add(remoteJunction);
                }
            }
        }
    }

    private void handleInternalDelete(final ReadGroup readGroup, final ReadRecord read, int minDeleteLength)
    {
        int maxDelete = maxDeleteLength(read.cigar());

        if(maxDelete < minDeleteLength)
            return;

        // convert the location of the internal delete into a junction
        int junctionStartPos = read.start() - 1;
        int junctionEndPos = 0;
        for(CigarElement element : read.cigar())
        {
            if(element.getOperator() == M)
            {
                junctionStartPos += element.getLength();
            }
            else if(element.getOperator() == D)
            {
                if(element.getLength() >= minDeleteLength)
                {
                    junctionEndPos = junctionStartPos + element.getLength() + 1;
                    break;
                }

                junctionStartPos += element.getLength();
            }
            else
            {
                continue;
            }
        }

        if(junctionEndPos <= junctionStartPos)
            return;

        JunctionData junctionStart = getOrCreateJunction(read, junctionStartPos, POS_ORIENT);
        junctionStart.JunctionGroups.add(readGroup);

        JunctionData junctionEnd = getOrCreateJunction(read, junctionEndPos, NEG_ORIENT);
        junctionEnd.JunctionGroups.add(readGroup);
    }

    private JunctionData getOrCreateJunction(final ReadRecord read, final byte orientation)
    {
        int junctionPosition = orientation == NEG_ORIENT ? read.start() : read.end();
        return getOrCreateJunction(read, junctionPosition, orientation);
    }

    private JunctionData getOrCreateJunction(final ReadRecord read, final int junctionPosition, final byte orientation)
    {
        // junctions are stored in ascending order to make finding them more efficient, especially for supporting reads
        int index = 0;

        // use a binary search if count exeeds some level, but with buckets < 1000 bases likely unnecessary
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.Position == junctionPosition)
            {
                if(junctionData.Orientation == orientation)
                    return junctionData;
            }
            else if(junctionData.Position > junctionPosition)
            {
                break;
            }

            ++index;
        }

        JunctionData junctionData = new JunctionData(junctionPosition, orientation, read);
        mJunctions.add(index, junctionData);
        return junctionData;
    }

    private void assignSupportingReads(final ReadFilterConfig filterConfig)
    {
        // attempt to assign reads which support a junction, otherwise discard to keep to move to next bucket
        mInitialSupportingReadCount = mSupportingReads.size();

        int index = 0;
        while(index < mSupportingReads.size())
        {
            ReadRecord read = mSupportingReads.get(index);

            if(!readSupportsJunction(read, filterConfig))
            {
                if(read.end() > mRegion.end())
                {
                    ++index; // will test support in the next bucket
                    continue;
                }
            }

            mSupportingReads.remove(index);
        }
    }

    private boolean readSupportsJunction(final ReadRecord read, final ReadFilterConfig filterConfig)
    {
        boolean supportsJunction = false;

        if(mJunctions.size() < 20)
        {
            for(JunctionData junctionData : mJunctions)
            {
                // any length soft clipping at the same base
                if(supportsJunction(read, junctionData, filterConfig))
                {
                    junctionData.SupportingReads.add(read);
                    supportsJunction = true;
                }
            }

            return supportsJunction;
        }

        int closeJunctionIndex = findJunctionIndex(read);

        if(closeJunctionIndex < 0)
            return false;

        // check up and down from this location
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);

            int index = searchUp ? closeJunctionIndex + 1 : closeJunctionIndex - 1;

            while(index >= 0 && index < mJunctions.size())
            {
                JunctionData junctionData = mJunctions.get(index);

                if(!readWithinJunctionRange(read, junctionData))
                    break;

                if(!reachedSupportingReadLimit(junctionData, filterConfig.MaxJunctionSupportingReads)) // to limit processing
                {
                    if(supportsJunction(read, junctionData, filterConfig))
                    {
                        junctionData.SupportingReads.add(read);
                        supportsJunction = true;
                    }
                }

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }

        return supportsJunction;
    }

    private static boolean reachedSupportingReadLimit(final JunctionData junctionData, int maxSupportingReads)
    {
        return maxSupportingReads > 0 && junctionData.supportingReadCount() >= maxSupportingReads;
    }

    private int findJunctionIndex(final ReadRecord read)
    {
        // binary search on junctions if the collection gets too large
        int currentIndex = mJunctions.size() / 2;
        int lowerIndex = 0;
        int upperIndex = mJunctions.size() - 1;

        while(true)
        {
            JunctionData junctionData = mJunctions.get(currentIndex);

            if(readWithinJunctionRange(read, junctionData))
                return currentIndex;

            if(positionWithin(read.start(), read.end(), junctionData.Position))
                return -1;

            if(read.end() < junctionData.Position)
            {
                // search lower
                if(lowerIndex + 1 == currentIndex)
                    return currentIndex;

                upperIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else if(read.start() > junctionData.Position)
            {
                // search higher
                if(currentIndex + 1 == upperIndex)
                    return currentIndex;

                lowerIndex = currentIndex;
                currentIndex = (lowerIndex + upperIndex) / 2;
            }
            else
            {
                break;
            }
        }

        return -1;
    }

    private static boolean readWithinJunctionRange(final ReadRecord read, final JunctionData junctionData)
    {
        boolean rightClipped = read.cigar().isRightClipped();
        boolean leftClipped = read.cigar().isLeftClipped();
        if(rightClipped || !leftClipped)
        {
            if(abs(read.end() - junctionData.Position) <= SUPPORTING_READ_DISTANCE)
                return true;
        }

        if(leftClipped || !rightClipped)
        {
            if(abs(read.start() - junctionData.Position) <= SUPPORTING_READ_DISTANCE)
                return true;
        }

        return false;
    }

    private static boolean hasDiscordantSupport(final ReadRecord read, final JunctionData junctionData, int maxDistance)
    {
        // must have both positions leading up to but not past the junction and one of its remote junctions
        if(junctionData.RemoteJunctions.isEmpty())
            return false;

        if(junctionData.Orientation != read.orientation())
            return false;

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(read.end() > junctionData.Position || abs(read.end() - junctionData.Position) > maxDistance)
                return false;
        }
        else
        {
            if(read.start() < junctionData.Position || abs(read.end() - junctionData.Position) > maxDistance)
                return false;
        }

        int otherPosition;
        String otherChromosome;
        byte otherOrientation;

        if(read.hasSuppAlignment())
        {
            otherChromosome = read.supplementaryAlignment().Chromosome;
            otherPosition = read.supplementaryAlignment().Position;
            otherOrientation = read.supplementaryAlignment().Strand == SUPP_POS_STRAND ? POS_ORIENT : NEG_ORIENT;
        }
        else
        {
            otherChromosome = read.MateChromosome;
            otherOrientation = read.mateOrientation();
            otherPosition = read.MatePosStart;
        }

        for(RemoteJunction remoteJunction : junctionData.RemoteJunctions)
        {
            if(!remoteJunction.Chromosome.equals(otherChromosome))
                continue;

            if(remoteJunction.Orientation != otherOrientation)
                continue;

            if(abs(remoteJunction.Position - otherPosition) > maxDistance)
                continue;

            if(remoteJunction.Orientation == POS_ORIENT && otherPosition <= remoteJunction.Position)
            {
                ++remoteJunction.Support;
                return true;
            }
            else if(remoteJunction.Orientation == NEG_ORIENT && otherPosition >= remoteJunction.Position)
            {
                ++remoteJunction.Support;
                return true;
            }
        }

        return false;
    }

    public static boolean supportsJunction(final ReadRecord read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        if(!read.cigar().isLeftClipped() && !read.cigar().isRightClipped())
            return hasDiscordantSupport(read, junctionData, filterConfig.MaxDiscordantFragmentDistance);

        /*
        eg if there is a candidate read with 101M50S at base 1000, how does a supporting read need to align and match?
         In this case the soft clip is at 1050. Check all reads which have a soft clip on the same side within +/- 50 bases.
        eg. you may find a read with soft clip at 1047 with 10S.
        Check if the 10 bases match the last 3 M and first 7S of the original read allowing for low base qual mismatches.
         If they do then that counts as an additional read even though it is much shorter soft clip and does not exactly match the base.
         */
        final ReadRecord juncRead = junctionData.InitialRead;

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(!read.cigar().isRightClipped())
                return false;

            int readRightPos = read.end();

            if(readRightPos == junctionData.Position)
                return true;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readRightPos - junctionData.Position) > SUPPORTING_READ_DISTANCE)
                return false;

            int scLength = read.cigar().getLastCigarElement().getLength();
            int readLength = read.readBases().length();
            int readEndPosIndex = readLength - scLength - 1;

            int juncReadLength = juncRead.readBases().length();
            int juncReadScLength = juncRead.cigar().getLastCigarElement().getLength();
            int juncReadEndPosIndex = juncReadLength - juncReadScLength - 1;
            int endPosDiff = juncRead.end() - readRightPos;

            int junctionReadOffset = juncReadEndPosIndex - readEndPosIndex - endPosDiff;

            // test against all the read's right soft-clipped bases
            for(int i = readLength - scLength; i < readLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                return false;
            }
        }
        else
        {
            if(!read.cigar().isLeftClipped())
                return false;

            int readLeftPos = read.start();

            if(readLeftPos == junctionData.Position)
                return true;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readLeftPos - junctionData.Position) > SUPPORTING_READ_DISTANCE)
                return false;

            // test for a base match for the read's soft-clipped bases, allow for low-qual matches

            // read: SC length -> start position
            // junc: SC length -> start position
            // junc read index = sc length diff - position diff

            int scLength = read.cigar().getFirstCigarElement().getLength();

            int juncReadScLength = juncRead.cigar().getFirstCigarElement().getLength();
            int posOffset = juncRead.start() - readLeftPos;
            int softClipDiff = juncReadScLength - scLength;
            int junctionReadOffset = softClipDiff - posOffset;
            int juncReadLength = juncRead.readBases().length();

            for(int i = 0; i < scLength; ++i)
            {
                char readBase = read.readBases().charAt(i);

                int juncIndex = i + junctionReadOffset;
                if(juncIndex < 0 || juncIndex >= juncReadLength)
                    return false;

                char juncReadBase = juncRead.readBases().charAt(juncIndex);

                if(readBase == juncReadBase)
                    continue;

                if(read.baseQualities()[i] < LOW_BASE_QUALITY || juncRead.baseQualities()[juncIndex] < LOW_BASE_QUALITY)
                    continue;

                return false;
            }
        }

        return true;
    }
}
