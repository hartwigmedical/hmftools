package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.AMP;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.DEL;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;

import org.jetbrains.annotations.NotNull;

public class RegionGeneEvents
{
    private static final EventRegionComparator sRegionComparator = new EventRegionComparator();
    private final ChrBaseRegion mRegion;
    private final GermlineStatus mStatus;
    private final List<String> mGenes;

    public RegionGeneEvents(GermlineAmpDel ampDel)
    {
        mRegion = getChrBaseRegion(ampDel);
        mStatus = ampDel.NormalStatus;
        mGenes = new ArrayList<>();
        mGenes.add(ampDel.GeneName);
    }

    public HumanChromosome chromosome()
    {
        return mRegion.humanChromosome();
    }

    public boolean offer(GermlineAmpDel ampDel)
    {
        ChrBaseRegion region = getChrBaseRegion(ampDel);
        int regionComparison = sRegionComparator.compare(region, mRegion);
        if(regionComparison < 0)
        {
            String msg = String.format("Region %s precedes this region %s", region, mRegion);
            throw new IllegalArgumentException(msg);
        }
        else if(regionComparison == 0)
        {
            int statusComparison = ampDel.NormalStatus.compareTo(mStatus);
            if(statusComparison < 0)
            {
                String msg = String.format("Status %s precedes this region status %s", ampDel.NormalStatus, mStatus);
                throw new IllegalArgumentException(msg);
            }
            else if(statusComparison == 0)
            {
                mGenes.add(ampDel.GeneName);
                return true;
            }
        }
        return false;
    }

    public List<String> genes()
    {
        return mGenes;
    }

    public ChrBaseRegion region()
    {
        return mRegion;
    }

    public AmpDelRegionFrequency.EventType eventType()
    {
        return mStatus == GermlineStatus.AMPLIFICATION ? AMP : DEL;
    }

    @NotNull
    private static ChrBaseRegion getChrBaseRegion(final GermlineAmpDel ampDel)
    {
        return new ChrBaseRegion(ampDel.Chromosome, ampDel.RegionStart, ampDel.RegionEnd);
    }
}
