package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.LOW_BASE_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.SUPPORTING_READ_DISTANCE;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.maxIndelLength;
import static com.hartwig.hmftools.svprep.reads.ReadType.JUNCTION;
import static com.hartwig.hmftools.svprep.reads.ReadType.NO_SUPPORT;
import static com.hartwig.hmftools.svprep.reads.ReadType.SUPPORT;
import static com.hartwig.hmftools.svprep.reads.RemoteJunction.addRemoteJunction;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SoftClipSide;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;
import com.hartwig.hmftools.svprep.HotspotCache;

import htsjdk.samtools.CigarElement;

public class JunctionTracker
{
    private final ChrBaseRegion mRegion;
    private final ReadFilterConfig mFilterConfig;
    private final HotspotCache mHotspotCache;
    private final List<BaseRegion> mBlacklistRegions;

    private final Map<String,ReadGroup> mReadGroups; // keyed by readId
    private final List<JunctionData> mJunctions; // ordered by position
    private final List<ReadGroup> mJunctionGroups; // groups used to form a junction
    private final List<ReadGroup> mSupportingGroups; // groups supporing a junction
    private int mInitialSupportingFrags;

    public JunctionTracker(
            final ChrBaseRegion region, final ReadFilterConfig config, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
    {
        mRegion = region;
        mFilterConfig = config;
        mHotspotCache = hotspotCache;

        mBlacklistRegions = Lists.newArrayList();

        List<BaseRegion> chrRegions = blacklist.getRegions(mRegion.Chromosome);

        if(chrRegions != null)
        {
            chrRegions.stream().filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                    .forEach(x -> mBlacklistRegions.add(x));
        }

        mReadGroups = Maps.newHashMap();
        mJunctions = Lists.newArrayList();
        mJunctionGroups = Lists.newArrayList();
        mSupportingGroups = Lists.newArrayList();
        mInitialSupportingFrags = 0;
    }

    public List<JunctionData> junctions() { return mJunctions; }
    public List<ReadGroup> junctionGroups() { return mJunctionGroups; }
    public List<ReadGroup> supportingGroups() { return mSupportingGroups; }
    public int initialSupportingFrags() { return mInitialSupportingFrags; }

    public void processRead(final ReadRecord read)
    {
        ReadGroup readGroup = mReadGroups.get(read.id());

        if(readGroup == null)
        {
            mReadGroups.put(read.id(), new ReadGroup(read));
            return;
        }

        readGroup.addRead(read);

        if(readGroup.isSimpleComplete() && readGroup.allNoSupport()) // purge irrelevant groups
            mReadGroups.remove(readGroup);
    }

    public void createJunctions()
    {
        List<ReadGroup> candidateSupportGroups = Lists.newArrayList();

        for(ReadGroup readGroup : mReadGroups.values())
        {
            if(overlapsBlacklist(readGroup)) // avoid creating any junction from a blacklist region
                continue;

            // ignore any group with a short overlapping fragment, likely adapter
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), INSERT_MAP_OVERLAP)))
                continue;

            if(readGroup.reads().stream().allMatch(x -> x.readType() == NO_SUPPORT)) // ignore groups with only fully-filtered reads
                continue;

            if(readGroup.isIncomplete())
            {
                // still assign but suggests something has been missed
                SV_LOGGER.debug("readGroup({}) incomplete post partition slicing", readGroup);
            }

            if(isJunctionFragment(readGroup))
            {
                createJunction(readGroup);
            }
            else
            {
                candidateSupportGroups.add(readGroup);
            }
        }

        mInitialSupportingFrags = candidateSupportGroups.size();

        Set<JunctionData> supportedJunctions = Sets.newHashSet();
        for(ReadGroup readGroup : candidateSupportGroups)
        {
            supportedJunctions.clear();

            for(ReadRecord read : readGroup.reads())
            {
                if(read.readType() != NO_SUPPORT)
                    checkJunctionSupport(read, supportedJunctions);
            }

            if(!supportedJunctions.isEmpty())
            {
                supportedJunctions.forEach(x -> x.SupportingGroups.add(readGroup));
                supportedJunctions.forEach(x -> readGroup.addJunctionPosition(x.Position));
                mSupportingGroups.add(readGroup);
            }
        }
    }

    private boolean isJunctionFragment(final ReadGroup readGroup)
    {
        return readGroup.reads().stream().anyMatch(x -> x.readType() == JUNCTION);
    }

    private void createJunction(final ReadGroup readGroup)
    {
        List<JunctionData> junctions = Lists.newArrayList();
        List<RemoteJunction> remoteJunctions = Lists.newArrayList();

        for(ReadRecord read : readGroup.reads())
        {
            if(read.hasSuppAlignment())
            {
                addRemoteJunction(remoteJunctions, RemoteJunction.fromSupplementaryData(read.supplementaryAlignment()));
            }

            handleInternalIndel(readGroup, read);

            SoftClipSide scSide = SoftClipSide.fromCigar(read.cigar());

            if(scSide == null || scSide.Length < mFilterConfig.MinSoftClipLength)
                continue;

            byte orientation = scSide.isLeft() ? NEG_ORIENT : POS_ORIENT;
            int position = scSide.isLeft() ? read.start() : read.end();

            if(!mRegion.containsPosition(position))
            {
                // will only be cached if no junction is within this region
                addRemoteJunction(remoteJunctions, new RemoteJunction(mRegion.Chromosome, position, orientation));
            }
            else
            {
                JunctionData junctionData = getOrCreateJunction(read, orientation);

                if(!junctions.contains(junctionData))
                    junctions.add(junctionData);
            }
        }

        if(junctions.isEmpty())
            return;

        mJunctionGroups.add(readGroup);
        junctions.forEach(x -> x.JunctionGroups.add(readGroup));
        junctions.forEach(x -> readGroup.addJunctionPosition(x.Position));

        if(!remoteJunctions.isEmpty())
        {
            for(RemoteJunction remoteJunction : remoteJunctions)
            {
                for(JunctionData junctionData : junctions)
                {
                    // ignore remotes (typically) supplementaries which point to junction in another read in this group
                    if(remoteJunction.matches(mRegion.Chromosome, junctionData.Position, junctionData.Orientation))
                        continue;

                    addRemoteJunction(junctionData.RemoteJunctions, remoteJunction);
                }
            }
        }
    }

    private boolean overlapsBlacklist(final ReadGroup readGroup)
    {
        for(BaseRegion region : mBlacklistRegions)
        {
            if(readGroup.reads().stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), region.start(), region.end())))
                return true;
        }

        return false;
    }

    private void handleInternalIndel(final ReadGroup readGroup, final ReadRecord read)
    {
        int maxDelete = maxIndelLength(read.cigar());

        if(maxDelete < mFilterConfig.MinIndelLength)
            return;

        // convert the location of the internal delete or insert into a junction
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
                if(element.getLength() >= mFilterConfig.MinIndelLength)
                {
                    junctionEndPos = junctionStartPos + element.getLength() + 1;
                    break;
                }

                junctionStartPos += element.getLength();
            }
            else if(element.getOperator() == I)
            {
                if(element.getLength() >= mFilterConfig.MinIndelLength)
                {
                    junctionEndPos = junctionStartPos + 1;
                    break;
                }
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

    private void checkJunctionSupport(final ReadRecord read, final Set<JunctionData> supportedJunctions)
    {
        if(mJunctions.size() < 20)
        {
            for(JunctionData junctionData : mJunctions)
            {
                if(supportedJunctions.contains(junctionData))
                    continue;

                // any length soft clipping at the same base
                if(supportsJunction(read, junctionData, mFilterConfig))
                {
                    read.setReadType(SUPPORT);
                    supportedJunctions.add(junctionData);
                }
            }

            return;
        }

        int closeJunctionIndex = findJunctionIndex(read);

        if(closeJunctionIndex < 0)
            return;

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

                if(!reachedSupportingReadLimit(junctionData, mFilterConfig.MaxJunctionSupportingReads)) // to limit processing
                {
                    if(!supportedJunctions.contains(junctionData) && supportsJunction(read, junctionData, mFilterConfig))
                    {
                        read.setReadType(SUPPORT);
                        supportedJunctions.add(junctionData);
                    }
                }

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }
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

    public void filterJunctions()
    {
        // alternatively when processing then also just remove all processed junctions
        int index = 0;
        while(index < mJunctions.size())
        {
            JunctionData junctionData = mJunctions.get(index);

            if(junctionData.totalSupport() >= mFilterConfig.MinJunctionSupport)
            {
                ++index;
                continue;
            }

            // check for a hotspot match
            boolean matchesHotspot = junctionData.RemoteJunctions.stream()
                    .anyMatch(x -> mHotspotCache.matchesHotspot(
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
}
