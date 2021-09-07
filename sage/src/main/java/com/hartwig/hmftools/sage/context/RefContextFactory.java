package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.count.EvictingArray;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.select.PanelSelector;

public class RefContextFactory
{
    private final SageConfig mConfig;
    private final String mSample;
    private final EvictingArray<RefContext> mRollingCandidates;
    private final PanelSelector mPanelSelector;
    private final List<AltContext> mSavedCandidates = Lists.newArrayList();

    public RefContextFactory(
            final SageConfig config, final String sample, final List<VariantHotspot> hotspots, final List<ChrBaseRegion> panel)
    {
        mSample = sample;
        mConfig = config;
        mPanelSelector = new PanelSelector(panel);

        final Predicate<AltContext> altContextPredicate = config.Filter.altContextFilter(new HotspotSelector(hotspots));

        final Consumer<RefContext> evictionHandler = (refContext) -> refContext.alts()
                .stream()
                .filter(AltContext::finaliseAndValidate)
                .filter(this::refPredicate)
                .filter(altContextPredicate)
                .forEach(mSavedCandidates::add);

        mRollingCandidates = new EvictingArray<>(256, evictionHandler);
    }

    public RefContext refContext(final String chromosome, final long position)
    {
        int maxDepth = maxReadDepth(chromosome, position);
        return mRollingCandidates.computeIfAbsent(position, aLong -> new RefContext(mSample, chromosome, position, maxDepth));
    }

    public List<AltContext> altContexts()
    {
        mRollingCandidates.evictAll();
        Collections.sort(mSavedCandidates);
        return mSavedCandidates;
    }

    private int maxReadDepth(final String chromosome, final long position)
    {
        return MitochondrialChromosome.contains(chromosome) || mPanelSelector.inPanel(position, position)
                ? mConfig.MaxReadDepthPanel : mConfig.MaxReadDepth;
    }

    private boolean refPredicate(final AltContext altContext)
    {
        for(int i = 0; i < altContext.ref().length(); i++)
        {
            char base = altContext.ref().charAt(i);
            if(base != 'G' && base != 'A' && base != 'T' && base != 'C')
            {
                return false;
            }
        }

        return true;
    }
}
