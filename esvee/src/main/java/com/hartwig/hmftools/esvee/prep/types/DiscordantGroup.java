package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class DiscordantGroup
{
    public final ChrBaseRegion Region;
    public final Orientation Orient;

    private PrepRead mInnermostRead;
    private final List<ReadGroup> mReadGroups;
    private final List<DiscordantRemoteRegion> mRemoteRegions;

    public DiscordantGroup(final ChrBaseRegion region, final Orientation orient, final ReadGroup readGroup, final PrepRead initialRead)
    {
        Region = region;
        Orient = orient;
        mReadGroups = Lists.newArrayList(readGroup);
        mInnermostRead = initialRead;
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

        Region.setEnd(max(Region.end(), read.end()));
        Region.setStart(min(Region.start(), read.start()));

        addRemoteRegionRead(readGroup, read);

        return true;
    }

    public void addRemoteRegionRead(final ReadGroup readGroup, final PrepRead read)
    {
        String remoteChromosome = read.MateChromosome;
        int remotePosStart = read.MatePosStart;
        int remotePosEnd = remotePosStart + read.record().getReadBases().length;
        Orientation remoteOrientation = read.mateOrientation();

        DiscordantRemoteRegion remoteRegion = mRemoteRegions.stream()
                .filter(x -> x.overlaps(remoteChromosome, remotePosStart, remotePosEnd))
                .filter(x -> x.Orient == remoteOrientation)
                .findFirst().orElse(null);

        if(remoteRegion != null)
        {
            remoteRegion.setStart(min(remoteRegion.start(), remotePosStart));
            remoteRegion.setEnd(max(remoteRegion.end(), remotePosEnd));
        }
        else
        {
            remoteRegion = new DiscordantRemoteRegion(
                    new ChrBaseRegion(remoteChromosome, remotePosStart, remotePosEnd), remoteOrientation);
            mRemoteRegions.add(remoteRegion);
        }

        remoteRegion.ReadGroups.add(readGroup);
    }

    public List<DiscordantRemoteRegion> remoteRegions() { return mRemoteRegions; }

    public static PrepRead firstPrimaryRead(final ReadGroup readGroup)
    {
        return readGroup.reads().stream().filter(x -> !x.record().getSupplementaryAlignmentFlag()).findFirst().orElse(null);
    }

    public static DiscordantGroup fromReadGroup(final ReadGroup readGroup)
    {
        PrepRead firstRead = firstPrimaryRead(readGroup);

        ChrBaseRegion region = new ChrBaseRegion(firstRead.Chromosome, firstRead.start(), firstRead.end());
        DiscordantGroup discordantGroup = new DiscordantGroup(region, firstRead.orientation(), readGroup, firstRead);

        discordantGroup.addRemoteRegionRead(readGroup, firstRead);

        return discordantGroup;
    }

    public String toString()
    {
        return format("region(%s) orient(%d) reads(%d) remoteRegions(%d)", Region, Orient.asByte(), mReadGroups.size(), mRemoteRegions.size());
    }
}
