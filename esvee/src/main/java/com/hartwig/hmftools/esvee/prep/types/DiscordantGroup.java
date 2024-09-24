package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class DiscordantGroup
{
    public final ChrBaseRegion Region;
    public final Orientation Orient;

    private PrepRead mInnermostRead;
    private final List<ReadGroup> mReadGroups;
    private final List<RemoteRegionReads> mRemoteRegions;

    public DiscordantGroup(final ChrBaseRegion region, final Orientation orient, final ReadGroup readGroup, final PrepRead initialRead)
    {
        Region = region;
        Orient = orient;
        mReadGroups = Lists.newArrayList(readGroup);
        mRemoteRegions = Lists.newArrayList();
    }

    public int innerPosition() { return Orient.isForward() ? Region.end() : Region.start(); }
    public PrepRead innerRead() { return mInnermostRead; }
    public List<ReadGroup> readGroups() { return mReadGroups; }

    public boolean tryAddReadGroup(final ReadGroup readGroup, final PrepRead read)
    {
        if(read.orientation() != Orient)
            return false;

        if(!Region.overlaps(read.Chromosome, read.start(), read.end()))
            return false;

        mReadGroups.add(readGroup);

        if(Orient.isForward())
            Region.setEnd(max(Region.end(), read.end()));
        else
            Region.setStart(min(Region.start(), read.start()));

        addRemoteRegionRead(read.MateChromosome, read.MatePosStart, read.MatePosStart + read.record().getReadBases().length);

        return true;
    }

    public void addRemoteRegionRead(final String chromosome, final int positionStart, final int positionEnd)
    {
        RemoteRegionReads existingRemote = mRemoteRegions.stream()
                .filter(x -> x.overlaps(chromosome, positionStart, positionEnd)).findFirst().orElse(null);

        if(existingRemote != null)
        {
            existingRemote.setStart(min(existingRemote.start(), positionStart));
            existingRemote.setEnd(max(existingRemote.end(), positionEnd));
            ++existingRemote.ReadCount;
        }
        else
        {
            mRemoteRegions.add(new RemoteRegionReads(new ChrBaseRegion(chromosome, positionStart, positionEnd)));
        }
    }

    public List<ChrBaseRegion> validRemoteRegions(int minReadCount)
    {
        return mRemoteRegions.stream().filter(x -> x.ReadCount >= minReadCount).collect(Collectors.toList());
    }

    private class RemoteRegionReads extends ChrBaseRegion
    {
        public int ReadCount;

        public RemoteRegionReads(final ChrBaseRegion region)
        {
            super(region.Chromosome, region.start(), region.end());
            ReadCount = 1;
        }

        public String toString() { return format("%s reads(%d)", super.toString(), ReadCount); }

    }

    public static PrepRead firstPrimaryRead(final ReadGroup readGroup)
    {
        return readGroup.reads().stream().filter(x -> !x.record().getSupplementaryAlignmentFlag()).findFirst().orElse(null);
    }

    public static DiscordantGroup fromReadGroup(final ReadGroup readGroup)
    {
        PrepRead firstRead = firstPrimaryRead(readGroup);

        ChrBaseRegion region = new ChrBaseRegion(firstRead.Chromosome, firstRead.start(), firstRead.end());
        DiscordantGroup discordantGroup = new DiscordantGroup(region, firstRead.orientation(), readGroup, firstRead);

        discordantGroup.addRemoteRegionRead(
                firstRead.MateChromosome, firstRead.MatePosStart, firstRead.MatePosStart + firstRead.record().getReadBases().length);

        return discordantGroup;
    }

    public String toString()
    {
        return format("region(%s) reads(%d) remoteRegions(%d)", Region, mReadGroups.size(), mRemoteRegions.size());
    }
}
