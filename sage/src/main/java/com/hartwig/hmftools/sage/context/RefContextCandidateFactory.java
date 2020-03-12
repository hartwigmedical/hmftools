package com.hartwig.hmftools.sage.context;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.count.EvictingArray;
import com.hartwig.hmftools.sage.select.HotspotSelector;

import org.jetbrains.annotations.NotNull;

public class RefContextCandidateFactory implements RefContextFactory {

    private final SageConfig config;
    private final HotspotSelector hotspotSelector;
    private final String sample;
    private final EvictingArray<RefContextCandidate> rollingCandidates;
    private final List<AltContext> savedCandidates = Lists.newArrayList();

    public RefContextCandidateFactory(@NotNull final SageConfig config, @NotNull final HotspotSelector hotspotSelector,
            @NotNull final String sample) {
        this.sample = sample;
        this.config = config;
        this.hotspotSelector = hotspotSelector;
        final Consumer<RefContextCandidate> evictionHandler = (refContext) -> {
            refContext.alts()
                    .stream()
                    .filter(this::refPredicate)
                    .filter(this::rawPredicate)
                    .filter(x -> x.primaryReadContext().readContext().isComplete())
                    .forEach(savedCandidates::add);
        };

        this.rollingCandidates = new EvictingArray<>(256, evictionHandler);
    }

    @NotNull
    public RefContextCandidate refContext(@NotNull final String chromosome, final long position) {
        return rollingCandidates.computeIfAbsent(position, aLong -> new RefContextCandidate(sample, chromosome, position));
    }

    @NotNull
    public List<AltContext> altContexts() {
        rollingCandidates.evictAll();
        Collections.sort(savedCandidates);
        return savedCandidates;
    }

    private boolean refPredicate(@NotNull final AltContext altContext) {
        for (int i = 0; i < altContext.ref().length(); i++) {
            char base = altContext.ref().charAt(i);
            if (base != 'G' && base != 'A' && base != 'T' && base != 'C') {
                return false;
            }
        }

        return true;
    }

    private boolean rawPredicate(@NotNull final AltContext altContext) {
        return hotspotSelector.isHotspot(altContext) || altContext.rawSupportAlt() >= config.filter().hardMinTumorRawAltSupport()
                && altContext.rawBaseQualityAlt() >= config.filter().hardMinTumorRawBaseQuality();
    }

}
