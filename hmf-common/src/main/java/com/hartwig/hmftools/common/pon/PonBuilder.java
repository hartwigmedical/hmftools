package com.hartwig.hmftools.common.pon;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PonBuilder {

    private final Map<VariantHotspot, Counter> map = Maps.newHashMap();

    public void add(@NotNull final VariantContext context) {
        final VariantHotspot hotspot = hotspot(context);
        final Counter counter = map.computeIfAbsent(hotspot, Counter::new);
        counter.increment();
    }

    @NotNull
    public List<VariantContext> build() {
        return map.values()
                .stream()
                .filter(x -> x.counter > 1)
                .sorted(Comparator.comparing(o -> o.hotspot))
                .map(PonBuilder::context)
                .collect(Collectors.toList());
    }

    @NotNull
    private static VariantHotspot hotspot(@NotNull final VariantContext context) {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(context.getAlternateAllele(0).getBaseString())
                .build();
    }

    @NotNull
    private static VariantContext context(@NotNull final Counter counter) {
        final Allele ref = Allele.create(counter.hotspot.ref(), true);
        final Allele alt = Allele.create(counter.hotspot.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        return new VariantContextBuilder().chr(counter.hotspot.chromosome())
                .start(counter.hotspot.position())
                .attribute(PonVCF.PON_COUNT, counter.counter)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, (int) counter.hotspot.position())
                .make();
    }

    static class Counter {
        private final VariantHotspot hotspot;
        private int counter;

        Counter(final VariantHotspot hotspot) {
            this.hotspot = hotspot;
        }

        void increment() {
            counter++;
        }
    }

}
