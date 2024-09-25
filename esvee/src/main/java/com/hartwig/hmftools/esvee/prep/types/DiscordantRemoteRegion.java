package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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
}
