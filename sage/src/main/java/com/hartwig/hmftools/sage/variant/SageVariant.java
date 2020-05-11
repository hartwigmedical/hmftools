package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class SageVariant implements GenomePosition {

    private final Set<String> filters;
    private final SageVariantTier tier;
    private final VariantHotspot variant;
    private final List<ReadContextCounter> normalAltContexts;
    private final List<ReadContextCounter> tumorAltContexts;

    private int localPhaseSet;
    private int localRealignSet;
    private int mixedImpact;
    private int phasedInframeIndel;
    private boolean realigned;

    public SageVariant(@NotNull final SageVariantTier tier, @NotNull final VariantHotspot variant, @NotNull final Set<String> filters,
            final List<ReadContextCounter> normal, final List<ReadContextCounter> tumorAltContexts) {
        assert (!tumorAltContexts.isEmpty());
        this.tier = tier;
        this.normalAltContexts = normal;
        this.tumorAltContexts = tumorAltContexts;
        this.filters = filters;
        this.variant = variant;
    }

    @NotNull
    public String ref() {
        return variant.ref();
    }

    @NotNull
    public String alt() {
        return variant.alt();
    }

    public int minNumberOfEvents() {
        return tumorAltContexts.get(0).minNumberOfEvents();
    }

    public boolean isRealigned() {
        return realigned;
    }

    public void realigned(final boolean realigned) {
        this.realigned = realigned;
    }

    public long end() {
        return position() + ref().length() - 1;
    }

    public boolean isIndel() {
        return variant.ref().length() != variant.alt().length();
    }

    public int phasedInframeIndel() {
        return phasedInframeIndel;
    }

    public void phasedInframeIndel(final int phasedInframeIndel) {
        this.phasedInframeIndel = phasedInframeIndel;
    }

    public boolean isMnv() {
        return variant.ref().length() >= 1 && variant.ref().length() == variant.alt().length();
    }

    public boolean isSnv() {
        return variant.ref().length() == 1 && variant.alt().length() == 1;
    }

    public boolean isInsert() {
        return variant.ref().length() < variant.alt().length();
    }

    public boolean isDelete() {
        return variant.ref().length() > variant.alt().length();
    }

    public int localPhaseSet() {
        return localPhaseSet;
    }

    public void localPhaseSet(int localPhaseSet) {
        this.localPhaseSet = localPhaseSet;
    }

    public int localRealignSet() {
        return localRealignSet;
    }

    public void localRealignSet(int localRealignSet) {
        this.localRealignSet = localRealignSet;
    }

    public int mixedGermlineImpact() {
        return mixedImpact;
    }

    public void mixedGermlineImpact(final int mixedImpact) {
        this.mixedImpact = mixedImpact;
    }

    public boolean isPassing() {
        return filters.isEmpty();
    }

    public boolean isTumorEmpty() {
        return tumorAltContexts.isEmpty();
    }

    public boolean isNormalEmpty() {
        return normalAltContexts.isEmpty();
    }

    @NotNull
    public VariantHotspot variant() {
        return variant;
    }

    @NotNull
    public SageVariantTier tier() {
        return tier;
    }

    @NotNull
    public Set<String> filters() {
        return filters;
    }

    @NotNull
    public ReadContext readContext() {
        return tumorAltContexts.get(0).readContext();
    }

    @NotNull
    public String microhomology() {
        return readContext().microhomology();
    }

    @NotNull
    public List<ReadContextCounter> normalAltContexts() {
        return normalAltContexts;
    }

    @NotNull
    public List<ReadContextCounter> tumorAltContexts() {
        return tumorAltContexts;
    }

    @NotNull
    @Override
    public String chromosome() {
        return variant().chromosome();
    }

    @Override
    public long position() {
        return variant().position();
    }

    public int totalQuality() {
        return tumorAltContexts.stream().mapToInt(ReadContextCounter::tumorQuality).sum();
    }

    public int maxQuality() {
        return tumorAltContexts.stream().mapToInt(ReadContextCounter::tumorQuality).max().orElse(0);
    }

}
