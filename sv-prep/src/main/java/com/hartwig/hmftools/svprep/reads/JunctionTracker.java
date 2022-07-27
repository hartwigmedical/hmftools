package com.hartwig.hmftools.svprep.reads;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvConstants.LOW_BASE_QUALITY;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;
import static com.hartwig.hmftools.svprep.SvConstants.MIN_MAP_QUALITY;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.POLY_G_SC;
import static com.hartwig.hmftools.svprep.reads.ReadFilters.isChimericRead;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.maxIndelLength;
import static com.hartwig.hmftools.svprep.reads.ReadType.EXACT_SUPPORT;
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
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.samtools.SoftClipSide;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;
import com.hartwig.hmftools.svprep.HotspotCache;
import com.hartwig.hmftools.svprep.SvConfig;

import htsjdk.samtools.CigarElement;

public class JunctionTracker
{
    private final ChrBaseRegion mRegion;
    private final ReadFilterConfig mFilterConfig;
    private final HotspotCache mHotspotCache;
    private final List<BaseRegion> mBlacklistRegions;

    private final Map<String,ReadGroup> mReadGroups; // keyed by readId
    private final Set<String> mExpectedReadIds; // as indicated by another partition

    private final List<JunctionData> mJunctions; // ordered by position
    private final List<ReadGroup> mJunctionGroups; // groups used to form a junction
    private final List<ReadGroup> mSupportingGroups; // groups supporting a junction
    private int mInitialSupportingFrags;
    private final int[] mBaseDepth;

    public JunctionTracker(
            final ChrBaseRegion region, final SvConfig config, final HotspotCache hotspotCache, final BlacklistLocations blacklist)
    {
        mRegion = region;
        mFilterConfig = config.ReadFiltering.config();
        mHotspotCache = hotspotCache;

        mBlacklistRegions = Lists.newArrayList();

        List<BaseRegion> chrRegions = blacklist.getRegions(mRegion.Chromosome);

        if(chrRegions != null)
        {
            chrRegions.stream().filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                    .forEach(x -> mBlacklistRegions.add(x));
        }

        mReadGroups = Maps.newHashMap();
        mExpectedReadIds = Sets.newHashSet();
        mJunctions = Lists.newArrayList();
        mJunctionGroups = Lists.newArrayList();
        mSupportingGroups = Lists.newArrayList();
        mInitialSupportingFrags = 0;

        mBaseDepth = config.CaptureDepth ? new int[mRegion.baseLength()] : null;
    }

    public List<JunctionData> junctions() { return mJunctions; }
    public List<ReadGroup> junctionGroups() { return mJunctionGroups; }
    public List<ReadGroup> supportingGroups() { return mSupportingGroups; }

    public List<ReadGroup> expectedGroups()
    {
        // groups not required by a junction but expected from other partitions
        return mReadGroups.values().stream()
                .filter(x -> x.isRemoteExpected() && x.junctionPositions() == null)
                .collect(Collectors.toList());
    }

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

        if(readGroup.isSimpleComplete()) // purge irrelevant groups
        {
            if(readGroup.allNoSupport())
            {
                captureDepth(readGroup);
                mReadGroups.remove(readGroup.id());
            }
            else if(groupInBlacklist(readGroup))
            {
                mReadGroups.remove(readGroup.id());
            }
        }
    }

    public void addExistingJunctions(final List<JunctionData> existingJunctions)
    {
        mJunctions.addAll(existingJunctions);
    }
    public void setExpectedReads(final Set<String> expectedReads) { mExpectedReadIds.addAll(expectedReads); }

    private void captureDepth(final ReadGroup readGroup)
    {
        if(mBaseDepth == null)
            return;

        int readsPosMin = readGroup.reads().stream().mapToInt(x -> x.start()).min().orElse(0);
        int readsPosMax = readGroup.reads().stream().mapToInt(x -> x.end()).max().orElse(0);
        int baseStart = max(readsPosMin - mRegion.start(), 0);
        int baseEnd = min(readsPosMax - mRegion.start(), mBaseDepth.length - 1);
        for(int i = baseStart; i <= baseEnd; ++i)
        {
            ++mBaseDepth[i];
        }
    }

    private int getBaseDepth(int position)
    {
        int baseIndex = position - mRegion.start();
        return baseIndex >= 0 && baseIndex < mBaseDepth.length ? mBaseDepth[baseIndex] : 0;
    }

    public void createJunctions()
    {
        List<ReadGroup> candidateSupportGroups = Lists.newArrayList();

        for(ReadGroup readGroup : mReadGroups.values())
        {
            captureDepth(readGroup);

            if(mExpectedReadIds.contains(readGroup.id()))
            {
                readGroup.markRemoteExpected();
                mExpectedReadIds.remove(readGroup.id());
            }

            // ignore any group with a short overlapping fragment, likely adapter
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), INSERT_MAP_OVERLAP)))
                continue;

            if(readGroup.reads().stream().allMatch(x -> x.readType() == NO_SUPPORT)) // ignore groups with only fully-filtered reads
                continue;

            // ignore any group with a poly-G insert (note that
            if(readGroup.reads().stream().anyMatch(x -> ReadFilterType.isSet(x.filters(), POLY_G_SC)))
                continue;

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

        if(mBaseDepth != null)
        {
            mJunctions.forEach(x -> x.setDepth(getBaseDepth(x.Position)));
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
            if(read.isUnmapped())
                continue;

            if(read.hasSuppAlignment())
            {
                addRemoteJunction(remoteJunctions, RemoteJunction.fromSupplementaryData(read.supplementaryAlignment()));
            }

            handleInternalIndel(readGroup, read);

            if(read.readType() != JUNCTION)
                continue;

            SoftClipSide scSide = SoftClipSide.fromCigar(read.cigar());

            if(scSide == null || scSide.Length < mFilterConfig.MinSoftClipLength)
                continue;

            byte orientation = scSide.isLeft() ? NEG_ORIENT : POS_ORIENT;
            int position = scSide.isLeft() ? read.start() : read.end();

            // junctions cannot fall in blacklist regions
            if(positionInBlacklist(position))
                continue;

            if(!mRegion.containsPosition(position))
            {
                addRemoteJunction(remoteJunctions, new RemoteJunction(mRegion.Chromosome, position, orientation));
            }
            else
            {
                JunctionData junctionData = getOrCreateJunction(read, orientation);

                if(!reachedFragmentCap(junctionData.junctionFragmentCount()) && !junctions.contains(junctionData))
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

                    junctionData.addRemoteJunction(remoteJunction);
                }
            }
        }
    }

    private boolean groupInBlacklist(final ReadGroup readGroup)
    {
        // test whether every read is in a blacklist region
        return readGroup.reads().stream().allMatch(x -> readInBlacklist(x));
    }

    private boolean readInBlacklist(final ReadRecord read)
    {
        return mBlacklistRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), read.start(), read.end()));
    }

    private boolean positionInBlacklist(int junctionPosition)
    {
        return mBlacklistRegions.stream().anyMatch(x -> x.containsPosition(junctionPosition));
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

        if(positionInBlacklist(junctionStartPos) || positionInBlacklist(junctionEndPos))
            return;

        JunctionData junctionStart = getOrCreateJunction(read, junctionStartPos, POS_ORIENT);
        JunctionData junctionEnd = getOrCreateJunction(read, junctionEndPos, NEG_ORIENT);

        junctionStart.markInternalIndel();
        junctionEnd.markInternalIndel();

        if(reachedFragmentCap(junctionStart.junctionFragmentCount()) && reachedFragmentCap(junctionEnd.junctionFragmentCount()))
            return;

        junctionStart.JunctionGroups.add(readGroup);
        junctionEnd.JunctionGroups.add(readGroup);

        mJunctionGroups.add(readGroup);
        readGroup.addJunctionPosition(junctionStartPos);
        readGroup.addJunctionPosition(junctionEndPos);
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
                checkJunctionSupport(read, junctionData, supportedJunctions);
            }

            return;
        }

        int closeJunctionIndex = findJunctionIndex(read);

        if(closeJunctionIndex < 0)
            return;

        checkJunctionSupport(read, mJunctions.get(closeJunctionIndex), supportedJunctions);

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

                checkJunctionSupport(read, junctionData, supportedJunctions);

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }
    }

    private void checkJunctionSupport(final ReadRecord read, final JunctionData junctionData, final Set<JunctionData> supportedJunctions)
    {
        if(reachedFragmentCap(junctionData.supportingFragmentCount())) // to limit processing
            return;

        if(supportedJunctions.contains(junctionData))
            return;

        // supporting reads cannot fall in blacklist regions
        if(readInBlacklist(read))
            return;

        if(hasExactJunctionSupport(read, junctionData, mFilterConfig))
        {
            read.setReadType(EXACT_SUPPORT);
            supportedJunctions.add(junctionData);
        }
        else if(hasDiscordantJunctionSupport(read, junctionData, mFilterConfig))
        {
            read.setReadType(SUPPORT);
            supportedJunctions.add(junctionData);
        }
    }

    private boolean reachedFragmentCap(final int fragments)
    {
        return mFilterConfig.JunctionFragmentCap > 0 && fragments >= mFilterConfig.JunctionFragmentCap;
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

    private boolean readWithinJunctionRange(final ReadRecord read, final JunctionData junctionData)
    {
        boolean rightClipped = read.cigar().isRightClipped();
        boolean leftClipped = read.cigar().isLeftClipped();
        if(rightClipped || !leftClipped)
        {
            if(abs(read.end() - junctionData.Position) <= mFilterConfig.MaxDiscordantFragmentDistance)
                return true;
        }

        if(leftClipped || !rightClipped)
        {
            if(abs(read.start() - junctionData.Position) <= mFilterConfig.MaxDiscordantFragmentDistance)
                return true;
        }

        return false;
    }

    public static boolean hasDiscordantJunctionSupport(
            final ReadRecord read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        // correct orientation
        if(junctionData.Orientation != read.orientation())
            return false;

        // correct side of the junction
        int junctionDistance = 0;

        if(junctionData.Orientation == POS_ORIENT)
        {
            if(read.end() > junctionData.Position)
                return false;

            junctionDistance = abs(read.end() - junctionData.Position);
        }
        else
        {
            if(read.start() < junctionData.Position || abs(read.end() - junctionData.Position) > filterConfig.MaxDiscordantFragmentDistance)
                return false;

            junctionDistance = abs(read.start() - junctionData.Position);
        }

        // any soft-clipping on the correct side if close to the junction
        if(junctionDistance <= filterConfig.MinSupportingReadDistance)
        {
            if(junctionData.Orientation == POS_ORIENT && read.record().getCigar().isRightClipped())
                return true;

            if(junctionData.Orientation == NEG_ORIENT && read.record().getCigar().isLeftClipped())
                return true;
        }

        // otherwise can be distant if chimeric
        if(junctionDistance <= filterConfig.MaxDiscordantFragmentDistance)
            return isChimericRead(read.record(), filterConfig);

        return false;
    }

    public static boolean hasExactJunctionSupport(
            final ReadRecord read, final JunctionData junctionData, final ReadFilterConfig filterConfig)
    {
        if(!read.cigar().isLeftClipped() && !read.cigar().isRightClipped())
            return false;

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

            if(juncRead == null)
                return false;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readRightPos - junctionData.Position) > filterConfig.MinSupportingReadDistance)
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

            if(juncRead == null)
                return false;

            // within 50 bases with exact sequence match in between the soft clip locations
            if(abs(readLeftPos - junctionData.Position) > filterConfig.MinSupportingReadDistance)
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

            if(junctionHasSupport(junctionData))
                ++index;
            else
                mJunctions.remove(index);
        }
    }

    private boolean junctionHasSupport(final JunctionData junctionData)
    {
        // first deal with junctions loaded from another sample - keep these if they've found any possible support
        if(junctionData.isExisting())
            return !junctionData.JunctionGroups.isEmpty() || !junctionData.SupportingGroups.isEmpty();

        // 1 junction read, 3 exact supporting reads altogether and 1 map-qual read
        int junctionFrags = junctionData.JunctionGroups.size();

        boolean hasPassingMapQualRead = junctionData.JunctionGroups.stream().anyMatch(x -> x.hasTypeAndMapQual(JUNCTION, MIN_MAP_QUALITY));

        if(hasPassingMapQualRead && junctionFrags >= mFilterConfig.MinJunctionSupport)
            return true;

        // look in the exact matches for additional support
        if(!hasPassingMapQualRead)
        {
            hasPassingMapQualRead = junctionData.SupportingGroups.stream().anyMatch(x -> x.hasTypeAndMapQual(EXACT_SUPPORT, MIN_MAP_QUALITY));
        }

        if(!hasPassingMapQualRead)
            return false;

        int exactSupportCount = (int) junctionData.SupportingGroups.stream()
                .filter(x -> x.reads().stream().anyMatch(y -> y.readType() == EXACT_SUPPORT)).count();

        if(junctionFrags + exactSupportCount >= mFilterConfig.MinJunctionSupport)
            return true;

        // check for a hotspot match
        boolean matchesHotspot = junctionData.RemoteJunctions.stream()
                .anyMatch(x -> mHotspotCache.matchesHotspot(
                        mRegion.Chromosome, x.Chromosome, junctionData.Position, x.Position, junctionData.Orientation, x.Orientation));

        if(matchesHotspot)
        {
            junctionData.markHotspot();

            if(junctionFrags + exactSupportCount >= MIN_HOTSPOT_JUNCTION_SUPPORT)
                return true;
        }

        return false;
    }
}
