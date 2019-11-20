package com.hartwig.hmftools.sage.variant;

import static com.hartwig.hmftools.sage.SageVCF.PASS;
import static com.hartwig.hmftools.sage.SageVCF.PHASE;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_DIFFERENCE;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_DISTANCE;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_MICRO_HOMOLOGY;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.SageVCF.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.sage.SageVCF.TIER;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public class SageVariantContextFactoryImpl {

    @NotNull
    public static VariantContext germlineOnly(@NotNull final SageVariant entry) {
        final AltContext normal = entry.normal();

        final List<Allele> alleles = createAlleles(normal);
        final Genotype normalGenotype = createGenotype(alleles, normal);
        final List<Genotype> genotypes = Collections.singletonList(normalGenotype);
        final ReadContextCounter normalCounter = normal.primaryReadContext();
        return createContext(entry, alleles, genotypes, normalCounter);
    }

    @NotNull
    public static VariantContext pairedTumorNormal(@NotNull final SageVariant entry) {
        final AltContext normal = entry.normal();
        final List<AltContext> tumorContexts = entry.tumorAltContexts();

        assert (tumorContexts.size() >= 1);

        final AltContext firstTumor = tumorContexts.get(0);

        final Allele ref = Allele.create(normal.ref(), true);
        final Allele alt = Allele.create(normal.alt(), false);
        final List<Allele> alleles = Lists.newArrayList(ref, alt);
        final Genotype normalGenotype = createGenotype(alleles, normal);

        final List<Genotype> genotypes = tumorContexts.stream().map(x -> createGenotype(alleles, x)).collect(Collectors.toList());
        genotypes.add(0, normalGenotype);

        final ReadContextCounter firstTumorCounter = firstTumor.primaryReadContext();
        return createContext(entry, alleles, genotypes, firstTumorCounter);
    }

    @NotNull
    private static VariantContext createContext(@NotNull final SageVariant variant, @NotNull final List<Allele> alleles,
            @NotNull final List<Genotype> genotypes, @NotNull final ReadContextCounter counter) {
        final VariantContextBuilder builder = new VariantContextBuilder().chr(variant.chromosome())
                .start(variant.position())
                .attribute(READ_CONTEXT, counter.toString())
                .attribute(READ_CONTEXT_DIFFERENCE, counter.readContext().distanceCigar())
                .attribute(READ_CONTEXT_DISTANCE, counter.readContext().distance())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, counter.readContextVaf())
                .log10PError(counter.quality() / -10d)
                .source("SAGE")
                .computeEndFromAlleles(alleles, (int) variant.position())
                .attribute(TIER, variant.tier())
                .genotypes(genotypes)
                .alleles(alleles)
                .filters(variant.filters());

        if (!counter.readContext().microhomology().isEmpty()) {
            builder.attribute(READ_CONTEXT_MICRO_HOMOLOGY, counter.readContext().microhomology());
        }

        if (counter.readContext().repeatCount() > 0) {
            builder.attribute(READ_CONTEXT_REPEAT_COUNT, counter.readContext().repeatCount())
                    .attribute(READ_CONTEXT_REPEAT_SEQUENCE, counter.readContext().repeat());
        }

        if (variant.localPhaseSet() > 0) {
            builder.attribute(PHASE, variant.localPhaseSet());
        }

        final VariantContext context = builder.make();
        if (context.isNotFiltered()) {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
    }

    @NotNull
    private static List<Allele> createAlleles(@NotNull final VariantHotspot variant) {
        final Allele ref = Allele.create(variant.ref(), true);
        final Allele alt = Allele.create(variant.alt(), false);
        return Lists.newArrayList(ref, alt);
    }

    @NotNull
    private static Genotype createGenotype(@NotNull final List<Allele> alleles, @NotNull final AltContext evidence) {
        ReadContextCounter readContextCounter = evidence.primaryReadContext();

        return new GenotypeBuilder(evidence.sample()).DP(evidence.readDepth())
                .AD(new int[] { evidence.refSupport(), evidence.altSupport() })
                .attribute(READ_CONTEXT_QUALITY, readContextCounter.qual())
                .attribute(READ_CONTEXT_COUNT, readContextCounter.rcc())
                .attribute(READ_CONTEXT_IMPROPER_PAIR, readContextCounter.improperPair())
                .alleles(alleles)
                .make();
    }

}
