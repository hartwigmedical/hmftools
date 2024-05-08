package com.hartwig.hmftools.sage.pon;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PonBuilder
{
    private static final int MIN_OUTPUT_COUNT = 2;
    private static final int MIN_INPUT_ALLELIC_DEPTH = 3;

    private final Map<VariantHotspot, Counter> mMap = new java.util.concurrent.ConcurrentHashMap<>();

    public void add(@NotNull final VariantContext context)
    {
        final VariantHotspot hotspot = hotspot(context);

        if(hotspot.ref().contains("N"))
            return;

        final Counter counter = mMap.computeIfAbsent(hotspot, Counter::new);
        final Genotype genotype = context.getGenotype(0);

        if(genotype.hasDP() && genotype.hasAD())
        {
            int[] adFields = genotype.getAD();
            int allelicDepth = adFields[1];
            if(allelicDepth >= MIN_INPUT_ALLELIC_DEPTH)
            {
                counter.increment(allelicDepth);
            }
        }
    }

    @NotNull
    public List<VariantContext> build()
    {
        return mMap.values()
                .stream()
                .filter(x -> x.counter() >= MIN_OUTPUT_COUNT)
                .sorted(Comparator.comparing(o -> o.hotspot))
                .map(PonBuilder::context)
                .collect(Collectors.toList());
    }

    @NotNull
    private static VariantHotspot hotspot(@NotNull final VariantContext context)
    {
        return ImmutableVariantHotspotImpl.builder()
                .chromosome(context.getContig())
                .position(context.getStart())
                .ref(context.getReference().getBaseString())
                .alt(context.getAlternateAllele(0).getBaseString())
                .build();
    }

    @NotNull
    private static VariantContext context(@NotNull final Counter counter)
    {
        final Allele ref = Allele.create(counter.hotspot.ref(), true);
        final Allele alt = Allele.create(counter.hotspot.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);

        return new VariantContextBuilder().chr(counter.hotspot.chromosome())
                .start(counter.hotspot.position())
                .attribute(PonVCF.PON_COUNT, counter.counter)
                .attribute(PonVCF.PON_TOTAL, counter.total)
                .attribute(PonVCF.PON_MAX, counter.max)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, (int) counter.hotspot.position())
                .make();
    }

    static class Counter
    {
        private final VariantHotspot hotspot;
        private final AtomicInteger counter = new AtomicInteger();
        private final AtomicInteger total = new AtomicInteger();
        private final AtomicInteger max = new AtomicInteger();

        Counter(final VariantHotspot hotspot)
        {
            this.hotspot = hotspot;
        }

        public int counter()
        {
            return counter.intValue();
        }

        void increment(int depth)
        {
            counter.incrementAndGet();
            total.addAndGet(depth);
            max.set(Integer.max(max.get(), depth));
        }
    }
}
