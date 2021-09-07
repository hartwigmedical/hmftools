package com.hartwig.hmftools.cobalt.ratio;

import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public class GCRatioSupplier
{
    private final GCMedianReadCount mTumorGCMedianReadCount;
    private final GCMedianReadCount mReferenceGCMedianReadCount;
    private final ListMultimap<Chromosome, ReadRatio> mTumorRatios;
    private final ListMultimap<Chromosome, ReadRatio> mReferenceRatios;

    public GCRatioSupplier(
            final Multimap<Chromosome, GCProfile> gcProfiles, final Multimap<Chromosome, CobaltCount> counts)
    {
        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.createImproved(gcProfiles);

        final GCRatioNormalization tumorRatiosBuilder = new GCRatioNormalization();
        final GCRatioNormalization referenceRatiosBuilder = new GCRatioNormalization();

        for(Chromosome chromosome : counts.keySet())
        {
            for(CobaltCount cobaltPosition : counts.get(chromosome))
            {
                final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(cobaltPosition);
                if(optionalGCProfile.isPresent())
                {
                    final GCProfile gcProfile = optionalGCProfile.get();
                    referenceRatiosBuilder.addPosition(chromosome, gcProfile, cobaltPosition.referenceReadCount());
                    tumorRatiosBuilder.addPosition(chromosome, gcProfile, cobaltPosition.tumorReadCount());
                }
            }
        }

        mReferenceGCMedianReadCount = referenceRatiosBuilder.gcMedianReadCount();
        mReferenceRatios = referenceRatiosBuilder.build(mReferenceGCMedianReadCount);

        mTumorGCMedianReadCount = tumorRatiosBuilder.gcMedianReadCount();
        mTumorRatios = tumorRatiosBuilder.build(mTumorGCMedianReadCount);
    }

    ListMultimap<Chromosome, ReadRatio> referenceRatios()
    {
        return mReferenceRatios;
    }

    GCMedianReadCount referenceGCMedianReadCount()
    {
        return mReferenceGCMedianReadCount;
    }

    ListMultimap<Chromosome, ReadRatio> tumorRatios()
    {
        return mTumorRatios;
    }

    GCMedianReadCount tumorGCMedianReadCount()
    {
        return mTumorGCMedianReadCount;
    }
}
