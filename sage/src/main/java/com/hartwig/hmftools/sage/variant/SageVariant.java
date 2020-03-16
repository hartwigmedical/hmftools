package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class SageVariant implements GenomePosition {

    private final Set<String> filters;
    private final SageVariantTier tier;
    private final List<ReadContextCounter> normalAltContexts;
    private final List<ReadContextCounter> tumorAltContexts;

    private int localPhaseSet;

    public SageVariant(final SageVariantTier tier, @NotNull final Set<String> filters, final List<ReadContextCounter> normal, final List<ReadContextCounter> tumorAltContexts) {
        assert (!tumorAltContexts.isEmpty());
        this.tier = tier;
        this.normalAltContexts = normal;
        this.tumorAltContexts = tumorAltContexts;
        this.filters = filters;
    }

    @NotNull
    public String ref() {
        return primaryNormal().alt();
    }

    @NotNull
    public String alt() {
        return primaryNormal().alt();
    }

    public long end() {
        return position() + ref().length() - 1;
    }

    public boolean isIndel() {
        return primaryNormal().ref().length() != primaryNormal().alt().length();
    }

    public boolean isInsert() {
        return primaryNormal().ref().length() < primaryNormal().alt().length();
    }

    public boolean isDelete() {
        return primaryNormal().ref().length() > primaryNormal().alt().length();
    }

    public int localPhaseSet() {
        return localPhaseSet;
    }

    public void localPhaseSet(int localPhaseSet) {
        this.localPhaseSet = localPhaseSet;
    }

    public boolean isPassing() {
        return filters.isEmpty();
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
    public ReadContextCounter primaryNormal() {
        return normalAltContexts.get(0);
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
    public ReadContextCounter primaryTumor() {
        return tumorAltContexts.get(0);
    }

    @NotNull
    @Override
    public String chromosome() {
        return primaryNormal().chromosome();
    }

    @Override
    public long position() {
        return primaryNormal().position();
    }

    public int totalQuality() {
        return tumorAltContexts.stream().mapToInt(ReadContextCounter::tumorQuality).sum();
    }

    public int maxQuality() {
        return tumorAltContexts.stream().mapToInt(ReadContextCounter::tumorQuality).max().orElse(0);
    }

}
