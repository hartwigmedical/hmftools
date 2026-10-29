package com.hartwig.hmftools.purple.tools;

import java.io.IOException;
import java.util.Comparator;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;

public class RegionEventExtractor
{
    private static final Comparator<GermlineAmpDel> GERMLINE_AMP_DEL_COMPARATOR =
            Comparator.comparing((GermlineAmpDel o) -> HumanChromosome.fromString(o.Chromosome))
                    .thenComparingInt(o -> o.RegionStart)
                    .thenComparingInt(o -> o.RegionEnd)
                    .thenComparing(o -> o.NormalStatus);
    private final ListMultimap<HumanChromosome, RegionGeneEvents> mEvents = ArrayListMultimap.create();

    public RegionEventExtractor(String ampDelFile) throws IOException
    {
        List<GermlineAmpDel> germlineAmpDels = GermlineAmpDel.read(ampDelFile);
        germlineAmpDels.sort(GERMLINE_AMP_DEL_COMPARATOR);
        RegionGeneEvents currentEvent = null;
        for(GermlineAmpDel event : germlineAmpDels)
        {
            if(currentEvent != null)
            {
                boolean merged = currentEvent.offer(event);
                if(!merged)
                {
                    mEvents.put(currentEvent.chromosome(), currentEvent);
                    currentEvent = new RegionGeneEvents(event);
                }
            }
            else
            {
                currentEvent = new RegionGeneEvents(event);
            }
        }
        if(currentEvent != null)
        {
            mEvents.put(currentEvent.chromosome(), currentEvent);
        }
    }

    public ListMultimap<HumanChromosome, RegionGeneEvents> events()
    {
        return mEvents;
    }
}
