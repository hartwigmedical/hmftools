package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.REMOTE_REGION_MERGE_MARGIN;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;

public class DiscordantRemoteRegion extends ChrBaseRegion
{
    public final Orientation Orient;
    public final List<ReadGroup> ReadGroups;

    public DiscordantRemoteRegion(final ChrBaseRegion region, final Orientation orient)
    {
        super(region.Chromosome, region.start(), region.end());
        Orient = orient;
        ReadGroups = Lists.newArrayList();
    }

    public int readCount() { return ReadGroups.size(); }

    public String toString()
    {
        return format("%s orient(%d) reads(%d)", super.toString(), Orient.asByte(), ReadGroups.size());
    }

    public static void mergeRemoteRegions(final List<DiscordantRemoteRegion> regions)
    {
        Collections.sort(regions, Comparator.comparing(x -> x));

        int distanceBuffer = REMOTE_REGION_MERGE_MARGIN;

        int index = 0;
        while(index < regions.size() - 1)
        {
            DiscordantRemoteRegion region = regions.get(index);

            int nextIndex = index + 1;
            while(nextIndex < regions.size())
            {
                DiscordantRemoteRegion nextRegion = regions.get(nextIndex);

                boolean mergeRegions = region.Chromosome.equals(nextRegion.Chromosome)
                        && region.Orient == nextRegion.Orient
                        && (positionsOverlap(region.start(), region.end(), nextRegion.start(), nextRegion.end())
                        || region.end() >= nextRegion.start() - distanceBuffer); // within close proximity


                if(mergeRegions)
                {
                    regions.remove(nextIndex);
                    region.setEnd(max(region.end(), nextRegion.end()));
                    region.ReadGroups.addAll(nextRegion.ReadGroups);
                }
                else
                {
                    ++nextIndex;
                }
            }

            ++index;
        }
    }
}
