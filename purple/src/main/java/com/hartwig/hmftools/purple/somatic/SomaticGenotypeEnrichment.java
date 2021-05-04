package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticGenotypeEnrichment implements VariantContextEnrichment
{
    enum SomaticGenotypeStatus
    {
        HET,
        HOM
    }

    @NotNull
    private static SomaticGenotypeStatus status(@NotNull AllelicDepth depth)
    {
        return depth.alleleReadCount() == depth.totalReadCount() ? SomaticGenotypeStatus.HOM : SomaticGenotypeStatus.HET;
    }

    private final String mGermlineSample;
    private final String mTumorSample;
    private final Consumer<VariantContext> mConsumer;

    public SomaticGenotypeEnrichment(final String germlineSample, final String tumorSample, final Consumer<VariantContext> consumer)
    {
        mGermlineSample = germlineSample;
        mTumorSample = tumorSample;
        mConsumer = consumer;
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        // AllelicDepth tumorDepth = AllelicDepth.fromGenotype(tumorGenotype);

        // Genotype germlineGenotype = context.getGenotype(mGermlineSample);
        // AllelicDepth germlineDepth = AllelicDepth.fromGenotype(germlineGenotype);
        // SomaticGenotypeStatus germlineStatus = status(germlineDepth);

        Allele refAllele = context.getReference();
        Allele altAllele = context.getAlternateAllele(0);

        // set the germline status
        Genotype germlineGT = context.getGenotype(mGermlineSample);

        List<Allele> germlineAlleles = Lists.newArrayList();
        germlineAlleles.add(refAllele);
        germlineAlleles.add(refAllele);

        Genotype germlineGenotype = new GenotypeBuilder(germlineGT).alleles(germlineAlleles).make();

        List<Genotype> updatedGenotypes = Lists.newArrayList();
        updatedGenotypes.add(germlineGenotype);

        // set the tumor status
        VariantContextDecorator variant = new VariantContextDecorator(context);

        Genotype tumorGT = context.getGenotype(mTumorSample);
        SomaticGenotypeStatus tumorStatus = variant.biallelic() ? SomaticGenotypeStatus.HOM : SomaticGenotypeStatus.HET;

        List<Allele> tumorAlleles = Lists.newArrayList();
        if(tumorStatus.equals(SomaticGenotypeStatus.HOM))
        {
            tumorAlleles.add(altAllele);
        }
        else
        {
            tumorAlleles.add(refAllele);
        }
        tumorAlleles.add(altAllele);

        final Genotype tumorGenotype = new GenotypeBuilder(tumorGT).alleles(tumorAlleles).make();
        updatedGenotypes.add(tumorGenotype);

        VariantContextBuilder builder = new VariantContextBuilder(context).genotypes(updatedGenotypes);
        mConsumer.accept(builder.make());
    }

    @Override
    public void flush() {}

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template)
    {
        return template;
    }
}
