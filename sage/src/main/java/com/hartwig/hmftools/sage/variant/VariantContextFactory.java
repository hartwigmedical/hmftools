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
import static com.hartwig.hmftools.sage.vcf.VariantVCF.STRAND_BIAS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.append.CandidateSerialization;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public final class VariantContextFactory
{
    private static final List<Allele> NO_CALL = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    public static VariantContext create(final SageVariant variant)
    {
        final List<Genotype> genotypes = Lists.newArrayList();
        for(int i = 0; i < variant.normalAltContexts().size(); i++)
        {
            ReadContextCounter normalContext = variant.normalAltContexts().get(i);
            genotypes.add(createGenotype(normalContext));
        }

        variant.tumorAltContexts().stream().map(VariantContextFactory::createGenotype).forEach(genotypes::add);
        return createContext(variant, genotypes);
    }

    private static VariantContext createContext(final SageVariant variant, final List<Genotype> genotypes)
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

        final VariantContext context = builder.make();
        if(context.isNotFiltered())
        {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
    }

    public static Genotype createGenotype(final ReadContextCounter counter)
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
                .attribute(STRAND_BIAS, counter.strandBias())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, counter.vaf())
                .alleles(NO_CALL)
                .make();
    }
}
