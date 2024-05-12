package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

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

    public void processVariant(final SomaticVariant variant)
    {
        VariantContext origContext = variant.context();
        Allele refAllele = origContext.getReference();

        if(origContext.getAlleles().size() < 2)
            return;

        Allele altAllele = origContext.getAlternateAllele(0);

        List<Genotype> updatedGenotypes = Lists.newArrayList();

        // set the germline status if present
        if(mGermlineSample != null && !mGermlineSample.isEmpty() && origContext.getGenotype(mGermlineSample) != null)
        {
            Genotype germlineGT = origContext.getGenotype(mGermlineSample);

            List<Allele> germlineAlleles = Lists.newArrayList();
            germlineAlleles.add(refAllele);
            germlineAlleles.add(refAllele);

            Genotype germlineGenotype = new GenotypeBuilder(germlineGT).alleles(germlineAlleles).make();

            updatedGenotypes.add(germlineGenotype);
        }

        // set the tumor status
        Genotype tumorGT = origContext.getGenotype(mTumorSample);
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

        VariantContextBuilder builder = new VariantContextBuilder(origContext).genotypes(updatedGenotypes);

        // remove any fields set by Pave (in regression testing only)
        if(variant.context().hasAttribute(REPORTED_FLAG))
            builder.rmAttribute(REPORTED_FLAG);

        VariantContext newContext = builder.make();
        variant.setContext(newContext);
    }
}
