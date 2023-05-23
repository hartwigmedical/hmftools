package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;

public class GcNormalizedRatioBuilder implements RatioBuilder
{
    private final GCMedianReadCount mGCMedianReadCount;
    private final ArrayListMultimap<Chromosome, ReadRatio> mGcRatios;

    public GcNormalizedRatioBuilder(final Multimap<Chromosome, GCProfile> gcProfiles,
            final Multimap<Chromosome, ReadRatio> counts,
            boolean useInterpolatedMedian)
    {
        CB_LOGGER.info("Applying ratio gc normalization");

        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.create(gcProfiles.values());

        final GCRatioNormalization gcRatioNormalization = new GCRatioNormalization(useInterpolatedMedian);

        for(Chromosome chromosome : counts.keySet())
        {
            for(ReadRatio readCount : counts.get(chromosome))
            {
                if (readCount.ratio() < 0.0)
                {
                    // skipped masked regions
                    continue;
                }

                final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(readCount);
                if(optionalGCProfile.isPresent())
                {
                    final GCProfile gcProfile = optionalGCProfile.get();
                    gcRatioNormalization.addPosition(chromosome, gcProfile, readCount.ratio());
                }
            }
        }

        mGCMedianReadCount = gcRatioNormalization.gcMedianReadCount();
        mGcRatios = gcRatioNormalization.build(mGCMedianReadCount);
    }

    @Override
    public ArrayListMultimap<Chromosome, ReadRatio> ratios()
    {
        return mGcRatios;
    }

    public GCMedianReadCount gcMedianReadCount()
    {
        return mGCMedianReadCount;
    }
}
