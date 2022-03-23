package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticGenotypeEnrichment
{
    enum SomaticGenotypeStatus
    {
        HET,
        HOM
    }

    private final String mGermlineSample;
    private final String mTumorSample;

    public SomaticGenotypeEnrichment(final String germlineSample, final String tumorSample)
    {
        mGermlineSample = germlineSample;
        mTumorSample = tumorSample;
    }

    public VariantContext processVariant(final VariantContext context)
    {
        Allele refAllele = context.getReference();

        if(context.getAlleles().size() < 2)
            return context;

        Allele altAllele = context.getAlternateAllele(0);

        List<Genotype> updatedGenotypes = Lists.newArrayList();

        // set the germline status if present
        if(mGermlineSample != null && !mGermlineSample.isEmpty() && context.getGenotype(mGermlineSample) != null)
        {
            Genotype germlineGT = context.getGenotype(mGermlineSample);

            List<Allele> germlineAlleles = Lists.newArrayList();
            germlineAlleles.add(refAllele);
            germlineAlleles.add(refAllele);

            Genotype germlineGenotype = new GenotypeBuilder(germlineGT).alleles(germlineAlleles).make();

            updatedGenotypes.add(germlineGenotype);
        }

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
        return builder.make();
    }
}
