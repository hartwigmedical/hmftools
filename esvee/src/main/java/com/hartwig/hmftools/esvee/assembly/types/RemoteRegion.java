package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REMOTE_REGION_WEAK_SUPP_PERCENT;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RemoteRegion extends ChrBaseRegion
{
    private final byte mOrientation;
    private final Set<String> mReadIds; // used to link with remote assemblies

    public static final int REMOTE_READ_TYPE_JUNCTION_MATE = 0;
    public static final int REMOTE_READ_TYPE_JUNCTION_SUPP = 1;
    public static final int REMOTE_READ_TYPE_DISCORDANT_READ = 2;

    private final int[] mReadTypeCount;
    private int mSoftClipMapQualTotal; // from reads with supplementaries

    public RemoteRegion(final ChrBaseRegion region, final byte orientation, final String readId, final int readType)
    {
        super(region.Chromosome, region.start(), region.end());
        mOrientation = orientation;
        mReadIds = Sets.newHashSet(readId);
        mReadTypeCount = new int[REMOTE_READ_TYPE_DISCORDANT_READ+1];
        ++mReadTypeCount[readType];
        mSoftClipMapQualTotal = 0;
    }

    public void addReadDetails(final String readId, final int posStart, final int posEnd, final int readType)
    {
        setStart(min(start(), posStart));
        setEnd(max(end(), posEnd));
        mReadIds.add(readId);
        ++mReadTypeCount[readType];
    }

    public byte orientation() { return mOrientation; }
    public boolean isForward() { return mOrientation == POS_ORIENT; }

    public Set<String> readIds() { return mReadIds; }
    public int readCount() { return mReadIds.size(); }

    public int[] readTypeCounts() { return mReadTypeCount; }

    public int nonSuppReadCount()
    {
        return mReadTypeCount[REMOTE_READ_TYPE_JUNCTION_MATE] + mReadTypeCount[REMOTE_READ_TYPE_DISCORDANT_READ];
    }

    public boolean isSuppOnlyRegion() { return nonSuppReadCount() == 0; }

    public int softClipMapQualTotal() { return mSoftClipMapQualTotal; }
    public void addSoftClipMapQual(int softClipLength, int mapQual) { mSoftClipMapQualTotal += softClipLength * mapQual; }

    public boolean matches(final RemoteRegion other) { return mOrientation == other.mOrientation && overlaps(other); }

    public boolean overlaps(final String otherChr, final int otherPosStart, final int otherPosEnd, final byte otherOrientation)
    {
        return mOrientation == otherOrientation && overlaps(otherChr, otherPosStart, otherPosEnd);
    }

    public String toString()
    {
        return format("%s orient(%d) reads(%d) counts(mate=%d supp=%d disc=%d) softClipMapQual(%d)",
                super.toString(), mOrientation, mReadIds.size(), mReadTypeCount[REMOTE_READ_TYPE_JUNCTION_MATE],
                mReadTypeCount[REMOTE_READ_TYPE_JUNCTION_SUPP], mReadTypeCount[REMOTE_READ_TYPE_DISCORDANT_READ], mSoftClipMapQualTotal);
    }

    public static void mergeRegions(final List<RemoteRegion> regions)
    {
        Collections.sort(regions, Comparator.comparing(x -> x));

        int index = 0;
        while(index < regions.size() - 1)
        {
            RemoteRegion region = regions.get(index);

            int nextIndex = index + 1;
            while(nextIndex < regions.size())
            {
                RemoteRegion nextRegion = regions.get(nextIndex);

                if(region.overlaps(nextRegion) && region.mOrientation == nextRegion.mOrientation)
                {
                    regions.remove(nextIndex);
                    region.setEnd(max(region.end(), nextRegion.end()));
                    region.readIds().addAll(nextRegion.readIds());

                    for(int i = 0; i < region.readTypeCounts().length; ++i)
                    {
                        region.readTypeCounts()[i] += nextRegion.readTypeCounts()[i];
                        region.addSoftClipMapQual(1, nextRegion.softClipMapQualTotal());
                    }
                }
                else
                {
                    ++nextIndex;
                }
            }

            ++index;
        }
    }

    public static void purgeWeakSupplementaryRegions(final List<RemoteRegion> regions)
    {
        int nonSuppSupport = 0;
        int maxSoftClipMapQual = 0;
        boolean hasSuppOnlyRegions = false;

        for(RemoteRegion region : regions)
        {
            hasSuppOnlyRegions |= region.isSuppOnlyRegion();
            nonSuppSupport += region.nonSuppReadCount();
            maxSoftClipMapQual = max(maxSoftClipMapQual, region.softClipMapQualTotal());
        }

        if(!hasSuppOnlyRegions)
            return;

        int index = 0;
        while(index < regions.size())
        {
            RemoteRegion region = regions.get(index);

            if(region.isSuppOnlyRegion())
            {
                if(region.softClipMapQualTotal() < REMOTE_REGION_WEAK_SUPP_PERCENT * maxSoftClipMapQual
                || region.readCount() < REMOTE_REGION_WEAK_SUPP_PERCENT * nonSuppSupport)
                {
                    regions.remove(index);
                    continue;
                }
            }

            ++index;
        }
    }
}
