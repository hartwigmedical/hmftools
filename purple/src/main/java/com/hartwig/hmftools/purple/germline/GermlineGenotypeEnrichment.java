package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import com.google.common.collect.Lists;
import org.apache.commons.math3.distribution.PoissonDistribution;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineGenotypeEnrichment
{
    public static final String LOW_VAF_FILTER = "LOW_VAF";

    enum GermlineGenotypeStatus
    {
        HET,
        HOM_ALT,
        LOW_VAF
    }

    private final String mGermlineSample;
    private final String mTumorSample;

    public GermlineGenotypeEnrichment(final String germlineSample, final String tumorSample)
    {
        mGermlineSample = germlineSample;
        mTumorSample = tumorSample;
    }

    public static GermlineGenotypeStatus status(AllelicDepth depth)
    {
        if(depth.AlleleReadCount == depth.TotalReadCount)
        {
            return GermlineGenotypeStatus.HOM_ALT;
        }

        boolean isHighVaf = Doubles.greaterThan(depth.AlleleReadCount, 0.75 * depth.TotalReadCount);
        if(isHighVaf)
        {
            return Doubles.lessThan(homPoisson(depth), 0.005) ? GermlineGenotypeStatus.HOM_ALT : GermlineGenotypeStatus.HET;
        }

        boolean isLowVaf = Doubles.lessThan(depth.AlleleReadCount, 0.3 * depth.TotalReadCount);
        if(isLowVaf)
        {
            return Doubles.lessThan(lowVafPoisson(depth), 0.002) ? GermlineGenotypeStatus.LOW_VAF : GermlineGenotypeStatus.HET;
        }

        return GermlineGenotypeStatus.HET;
    }

    public void processVariant(final GermlineVariant variant)
    {
        VariantContext context = variant.context();
        Genotype germlineGenotype = context.getGenotype(mGermlineSample);
        Genotype tumorGenotype = context.getGenotype(mTumorSample);
        AllelicDepth germlineDepth = AllelicDepth.fromGenotype(germlineGenotype);
        AllelicDepth tumorDepth = AllelicDepth.fromGenotype(tumorGenotype);

        GermlineGenotypeStatus germlineStatus = status(germlineDepth);
        GermlineGenotypeStatus tumorStatus = status(tumorDepth);
        GermlineGenotypeStatus status = combined(germlineStatus, tumorStatus);

        Allele refAllele = context.getReference();
        Allele altAllele = context.getAlternateAllele(0);
        List<Allele> genotypeAlleles = Lists.newArrayList();
        if(status.equals(GermlineGenotypeStatus.HOM_ALT))
        {
            genotypeAlleles.add(altAllele);
        }
        else
        {
            genotypeAlleles.add(refAllele);
        }
        genotypeAlleles.add(altAllele);

        final Genotype updatedGenotype = new GenotypeBuilder(germlineGenotype).alleles(genotypeAlleles).make();
        final List<Genotype> updatedGenotypes = Lists.newArrayList();
        context.getGenotypes().stream().filter(x -> !x.getSampleName().equals(mGermlineSample)).forEach(updatedGenotypes::add);
        updatedGenotypes.add(updatedGenotype);

        VariantContextBuilder builder = new VariantContextBuilder(context).genotypes(updatedGenotypes);
        if(status.equals(GermlineGenotypeStatus.LOW_VAF))
        {
            builder.filter(LOW_VAF_FILTER);
            variant.filters().remove(PASS);
            variant.filters().add(LOW_VAF_FILTER);
        }

        VariantContext newContext = builder.make();
        variant.setContext(newContext);
    }

    static GermlineGenotypeStatus combined(GermlineGenotypeStatus germline, GermlineGenotypeStatus tumor)
    {
        if(germline == GermlineGenotypeStatus.HOM_ALT)
        {
            return tumor == GermlineGenotypeStatus.HOM_ALT ? GermlineGenotypeStatus.HOM_ALT : GermlineGenotypeStatus.HET;
        }

        return germline;
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFFilterHeaderLine(LOW_VAF_FILTER, "Germline variant has very low allelic frequency"));
        return template;
    }

    private static double homPoisson(AllelicDepth depth)
    {
        return new PoissonDistribution(depth.TotalReadCount / 2d).cumulativeProbability(depth.TotalReadCount - depth.AlleleReadCount);
    }

    private static double lowVafPoisson(AllelicDepth depth)
    {
        return new PoissonDistribution(depth.TotalReadCount / 2d).cumulativeProbability(depth.AlleleReadCount);
    }
}
