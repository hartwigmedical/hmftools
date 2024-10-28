package com.hartwig.hmftools.esvee.assembly.types;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REMOTE_REGION_MERGE_MARGIN;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REMOTE_REGION_WEAK_SUPP_PERCENT;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.DISCORDANT;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteReadType.SUPPLEMENTARY;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class RemoteRegion extends ChrBaseRegion
{
    private final Set<String> mReadIds; // used to link with remote assemblies, and note is trimmed to match ID in SupportRead

    private final int[] mReadTypeCount;
    private int mSoftClipMapQualTotal; // from reads with supplementaries

    public RemoteRegion(final ChrBaseRegion region, final String readId, final RemoteReadType readType)
    {
        super(region.Chromosome, region.start(), region.end());

        mReadIds = Sets.newHashSet(readId);

        mReadTypeCount = new int[RemoteReadType.values().length];
        ++mReadTypeCount[readType.ordinal()];

        mSoftClipMapQualTotal = 0;
    }

    public void addReadDetails(final String readId, final int posStart, final int posEnd, final RemoteReadType readType)
    {
        setStart(min(start(), posStart));
        setEnd(max(end(), posEnd));

        mReadIds.add(readId);
        ++mReadTypeCount[readType.ordinal()];
    }

    public Set<String> readIds() { return mReadIds; }
    public int readCount() { return mReadIds.size(); }

    public int[] readTypeCounts() { return mReadTypeCount; }

    public int nonSuppReadCount()
    {
        return mReadTypeCount[JUNCTION_MATE.ordinal()] + mReadTypeCount[DISCORDANT.ordinal()];
    }

    public boolean isSuppOnlyRegion() { return nonSuppReadCount() == 0; }

    public int softClipMapQualTotal() { return mSoftClipMapQualTotal; }
    public void addSoftClipMapQual(int softClipLength, int mapQual) { mSoftClipMapQualTotal += softClipLength * mapQual; }

    public boolean matches(final RemoteRegion other) { return overlaps(other); }

    public boolean overlaps(final String otherChr, final int otherPosStart, final int otherPosEnd)
    {
        return super.overlaps(otherChr, otherPosStart, otherPosEnd);
    }

    public boolean overlapsAssembly(final JunctionAssembly assembly)
    {
        return super.overlaps(assembly.junction().Chromosome, assembly.minAlignedPosition(), assembly.maxAlignedPosition());
    }

    public boolean hasReadId(final String fullReadId)
    {
        return mReadIds.stream().anyMatch(x -> fullReadId.equals(x));
    }

    public String toString()
    {
        return format("%s reads(%d) counts(mate=%d supp=%d disc=%d) softClipMapQual(%d)",
                super.toString(), mReadIds.size(), mReadTypeCount[JUNCTION_MATE.ordinal()],
                mReadTypeCount[SUPPLEMENTARY.ordinal()], mReadTypeCount[DISCORDANT.ordinal()], mSoftClipMapQualTotal);
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

                boolean mergeRegions = region.Chromosome.equals(nextRegion.Chromosome)
                        && (positionsOverlap(region.start(), region.end(), nextRegion.start(), nextRegion.end())
                        || region.end() >= nextRegion.start() - REMOTE_REGION_MERGE_MARGIN); // within close proximity


                if(mergeRegions)
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
        // purge any supplementary-only remote region if it's soft-clip read map-qual < 10% of the top other regions same value, OR
        // if the number of reads it has is < 10% of the reads in a non-supp-only region
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
