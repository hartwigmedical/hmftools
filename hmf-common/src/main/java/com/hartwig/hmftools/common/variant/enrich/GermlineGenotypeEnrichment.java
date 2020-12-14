package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;

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

public class GermlineGenotypeEnrichment implements VariantContextEnrichment {

    public static final String LOW_VAF_FILTER = "LOW_VAF";

    enum GermlineGenotypeStatus {
        HET,
        HOM,
        LOW_VAF
    }

    @NotNull
    static GermlineGenotypeStatus status(@NotNull AllelicDepth depth) {
        if (depth.alleleReadCount() == depth.totalReadCount()) {
            return GermlineGenotypeStatus.HOM;
        }

        boolean isHighVaf = Doubles.greaterThan(depth.alleleReadCount(), 0.75 * depth.totalReadCount());
        if (isHighVaf) {
            return Doubles.lessThan(homPoisson(depth), 0.005) ? GermlineGenotypeStatus.HOM : GermlineGenotypeStatus.HET;
        }

        boolean isLowVaf = Doubles.lessThan(depth.alleleReadCount(), 0.3 * depth.totalReadCount());
        if (isLowVaf) {
            return Doubles.lessThan(lowVafPoisson(depth), 0.002) ? GermlineGenotypeStatus.LOW_VAF : GermlineGenotypeStatus.HET;
        }

        return GermlineGenotypeStatus.HET;
    }

    @NotNull
    private final String germlineSample;
    @NotNull
    private final Consumer<VariantContext> consumer;

    public GermlineGenotypeEnrichment(@NotNull final String germlineSample, @NotNull final Consumer<VariantContext> consumer) {
        this.germlineSample = germlineSample;
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        Genotype genotype = context.getGenotype(germlineSample);
        AllelicDepth depth = AllelicDepth.fromGenotype(genotype);

        GermlineGenotypeStatus status = status(depth);
        if (status.equals(GermlineGenotypeStatus.LOW_VAF)) {
            VariantContext variantContext = new VariantContextBuilder(context).filter(LOW_VAF_FILTER).make();
            consumer.accept(variantContext);
        } else {
            Allele refAllele = context.getReference();
            Allele altAllele = context.getAlternateAllele(0);
            List<Allele> genotypeAlleles = Lists.newArrayList();
            if (status.equals(GermlineGenotypeStatus.HET)) {
                genotypeAlleles.add(refAllele);
            } else {
                genotypeAlleles.add(altAllele);
            }
            genotypeAlleles.add(altAllele);

            final Genotype updatedGenotype = new GenotypeBuilder(genotype).alleles(genotypeAlleles).make();
            final List<Genotype> updatedGenotypes = Lists.newArrayList();
            context.getGenotypes().stream().filter(x -> !x.getSampleName().equals(germlineSample)).forEach(updatedGenotypes::add);
            updatedGenotypes.add(updatedGenotype);

            VariantContext variantContext = new VariantContextBuilder(context).genotypes(updatedGenotypes).make();
            consumer.accept(variantContext);
        }
    }

    @Override
    public void flush() {

    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFFilterHeaderLine(LOW_VAF_FILTER, "Germline variant has very low allelic frequency"));
        return template;
    }

    static double homPoisson(@NotNull AllelicDepth depth) {
        return new PoissonDistribution(depth.totalReadCount() / 2d).cumulativeProbability(depth.totalReadCount() - depth.alleleReadCount());
    }

    static double lowVafPoisson(@NotNull AllelicDepth depth) {
        return new PoissonDistribution(depth.totalReadCount() / 2d).cumulativeProbability(depth.alleleReadCount());
    }
}
