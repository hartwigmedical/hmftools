package com.hartwig.hmftools.sage.variant;

import static com.hartwig.hmftools.common.sage.SageMetaData.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.sage.SageMetaData.LOCAL_REALIGN_SET;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.PASS;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_ALLELIC_BASE_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_ALLELIC_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RIGHT_ALIGNED_MICROHOMOLOGY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.candidate.CandidateSerialization;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public final class SageVariantContextFactory
{
    private static final List<Allele> NO_CALL = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    @NotNull
    public static VariantContext addGenotype(@NotNull final VariantContext parent, @NotNull final List<ReadContextCounter> counters)
    {
        final VariantContextBuilder builder = new VariantContextBuilder(parent);
        final List<Genotype> genotypes = Lists.newArrayList(parent.getGenotypes());

        for(ReadContextCounter counter : counters)
        {
            Genotype genotype = createGenotype(counter);
            genotypes.add(genotype);
        }
        return builder.genotypes(genotypes).make();
    }

    @NotNull
    public static VariantContext create(@NotNull final SageVariant entry)
    {
        final List<Genotype> genotypes = Lists.newArrayList();
        for(int i = 0; i < entry.normalAltContexts().size(); i++)
        {
            ReadContextCounter normalContext = entry.normalAltContexts().get(i);
            genotypes.add(createGenotype(normalContext));
        }

        entry.tumorAltContexts().stream().map(SageVariantContextFactory::createGenotype).forEach(genotypes::add);
        return createContext(entry, genotypes);
    }

    @NotNull
    private static VariantContext createContext(@NotNull final SageVariant variant, @NotNull final List<Genotype> genotypes)
    {
        final VariantContextBuilder builder = CandidateSerialization.toContext(variant.candidate())
                .log10PError(variant.totalQuality() / -10d)
                .genotypes(genotypes)
                .filters(variant.filters());

        if(variant.localPhaseSet() > 0)
        {
            builder.attribute(LOCAL_PHASE_SET, variant.localPhaseSet());
        }

        if(variant.localRealignSet() > 0)
        {
            builder.attribute(LOCAL_REALIGN_SET, variant.localRealignSet());
        }

        if(variant.mixedGermlineImpact() > 0)
        {
            builder.attribute(MIXED_SOMATIC_GERMLINE, variant.mixedGermlineImpact());
        }

        if(variant.isRealigned())
        {
            builder.attribute(RIGHT_ALIGNED_MICROHOMOLOGY, true);
        }

        final VariantContext context = builder.make();
        if(context.isNotFiltered())
        {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
    }

    @NotNull
    private static Genotype createGenotype(@NotNull final ReadContextCounter counter)
    {
        return new GenotypeBuilder(counter.Sample).DP(counter.depth())
                .AD(new int[] { counter.refSupport(), counter.altSupport() })
                .attribute(READ_CONTEXT_QUALITY, counter.quality())
                .attribute(READ_CONTEXT_COUNT, counter.counts())
                .attribute(READ_CONTEXT_IMPROPER_PAIR, counter.improperPair())
                .attribute(READ_CONTEXT_JITTER, counter.jitter())
                .attribute(RAW_ALLELIC_DEPTH, new int[] { counter.rawRefSupport(), counter.rawAltSupport() })
                .attribute(RAW_ALLELIC_BASE_QUALITY, new int[] { counter.rawRefBaseQuality(), counter.rawAltBaseQuality() })
                .attribute(RAW_DEPTH, counter.rawDepth())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, counter.vaf())
                .alleles(NO_CALL)
                .make();
    }
}
