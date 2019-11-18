package com.hartwig.hmftools.sage;

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

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public class SageVariantContextFactory {

    @NotNull
    public static VariantContext create(@NotNull final SageVariant entry) {
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

        final VariantContextBuilder builder = new VariantContextBuilder().chr(normal.chromosome())
                .start(normal.position())
                .attribute(READ_CONTEXT, firstTumorCounter.toString())
                .attribute(READ_CONTEXT_DIFFERENCE, firstTumorCounter.readContext().distanceCigar())
                .attribute(READ_CONTEXT_DISTANCE, firstTumorCounter.readContext().distance())
                .attribute(TIER, entry.tier())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, firstTumorCounter.readContextVaf())
                .computeEndFromAlleles(alleles, (int) normal.position())
                .source("SAGE")
                .genotypes(genotypes)
                .alleles(alleles)
                .filters(entry.filters());

        if (!firstTumor.primaryReadContext().readContext().microhomology().isEmpty()) {
            builder.attribute(READ_CONTEXT_MICRO_HOMOLOGY, firstTumorCounter.readContext().microhomology());
        }

        if (firstTumor.primaryReadContext().readContext().repeatCount() > 0) {
            builder.attribute(READ_CONTEXT_REPEAT_COUNT, firstTumorCounter.readContext().repeatCount())
                    .attribute(READ_CONTEXT_REPEAT_SEQUENCE, firstTumorCounter.readContext().repeat());
        }

        if (entry.localPhaseSet() > 0) {
            builder.attribute(PHASE, entry.localPhaseSet());
        }

        final VariantContext context = builder.make();
        context.getCommonInfo().setLog10PError(firstTumorCounter.quality() / -10d);
        if (context.isNotFiltered()) {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
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
