package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.count.EvictingArray;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.select.PanelSelector;

import org.jetbrains.annotations.NotNull;

public class RefContextFactory
{

    private final SageConfig config;
    private final String sample;
    private final EvictingArray<RefContext> rollingCandidates;
    private final PanelSelector<GenomeRegion> panelSelector;
    private final List<AltContext> savedCandidates = Lists.newArrayList();

    public RefContextFactory(@NotNull final SageConfig config, @NotNull final String sample, final List<VariantHotspot> hotspots,
            final List<GenomeRegion> panel)
    {
        this.sample = sample;
        this.config = config;
        this.panelSelector = new PanelSelector<>(panel);
        final Predicate<AltContext> altContextPredicate = config.filter().altContextFilter(new HotspotSelector(hotspots));
        final Consumer<RefContext> evictionHandler = (refContext) -> refContext.alts()
                .stream()
                .filter(AltContext::finaliseAndValidate)
                .filter(this::refPredicate)
                .filter(altContextPredicate)
                .forEach(savedCandidates::add);

        this.rollingCandidates = new EvictingArray<>(256, evictionHandler);
    }

    @NotNull
    public RefContext refContext(@NotNull final String chromosome, final long position)
    {
        int maxDepth = maxReadDepth(chromosome, position);
        return rollingCandidates.computeIfAbsent(position, aLong -> new RefContext(sample, chromosome, position, maxDepth));
    }

    @NotNull
    public List<AltContext> altContexts()
    {
        rollingCandidates.evictAll();
        Collections.sort(savedCandidates);
        return savedCandidates;
    }

    private int maxReadDepth(final String chromosome, final long position)
    {
        return MitochondrialChromosome.contains(chromosome) || panelSelector.inPanel(position, position)
                ? config.maxReadDepthPanel()
                : config.maxReadDepth();
    }

    private boolean refPredicate(@NotNull final AltContext altContext)
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
