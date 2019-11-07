package com.hartwig.hmftools.sage.variant;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.sage.context.AltContext;

import org.jetbrains.annotations.NotNull;

public class SageVariant {

    private final AltContext normal;
    private final Set<String> filters;
    private final SageVariantTier tier;
    private final List<AltContext> tumorAltContexts;

    public SageVariant(final SageVariantTier tier, @NotNull final Set<String> filters, final AltContext normal,
            final List<AltContext> tumorAltContexts) {
        assert (!tumorAltContexts.isEmpty());
        this.tier = tier;
        this.normal = normal;
        this.tumorAltContexts = tumorAltContexts;
        this.filters = filters;
    }

    @Deprecated
    public SageVariant(final AltContext normal, final List<AltContext> tumorAltContexts) {
        assert (!tumorAltContexts.isEmpty());

        this.tier = SageVariantTier.WIDE;
        this.normal = normal;
        this.tumorAltContexts = tumorAltContexts;
        this.filters = Sets.newHashSet();
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
    public AltContext normal() {
        return normal;
    }

    @NotNull
    public List<AltContext> tumorAltContexts() {
        return tumorAltContexts;
    }

    @NotNull
    public AltContext primaryTumor() {
        return tumorAltContexts.get(0);
    }

}
