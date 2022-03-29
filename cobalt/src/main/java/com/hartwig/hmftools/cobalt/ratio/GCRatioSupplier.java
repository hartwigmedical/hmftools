package com.hartwig.hmftools.cobalt.ratio;

import java.util.Optional;
import java.util.function.Function;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;

public class GCRatioSupplier
{
    private final GCMedianReadCount mGCMedianReadCount;
    private final ListMultimap<Chromosome, ReadRatio> mGcRatios;

    public GCRatioSupplier(
            final Multimap<Chromosome, GCProfile> gcProfiles, final Multimap<Chromosome, CobaltCount> counts,
            Function<CobaltCount, Integer> readCountGetter)
    {
        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.createImproved(gcProfiles);

        final GCRatioNormalization gcRatioNormalization = new GCRatioNormalization();

        for(Chromosome chromosome : counts.keySet())
        {
            for(CobaltCount cobaltPosition : counts.get(chromosome))
            {
                final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(cobaltPosition);
                if(optionalGCProfile.isPresent())
                {
                    final GCProfile gcProfile = optionalGCProfile.get();
                    gcRatioNormalization.addPosition(chromosome, gcProfile, readCountGetter.apply(cobaltPosition));
                }
            }
        }

        mGCMedianReadCount = gcRatioNormalization.gcMedianReadCount();
        mGcRatios = gcRatioNormalization.build(mGCMedianReadCount);
    }

    ListMultimap<Chromosome, ReadRatio> gcRatios()
    {
        return mGcRatios;
    }

    GCMedianReadCount gcMedianReadCount()
    {
        return mGCMedianReadCount;
    }
}
