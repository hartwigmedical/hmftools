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
        HOM_ALT,
        LOW_VAF
    }

    @NotNull
    static GermlineGenotypeStatus status(@NotNull AllelicDepth depth) {
        if (depth.alleleReadCount() == depth.totalReadCount()) {
            return GermlineGenotypeStatus.HOM_ALT;
        }

        boolean isHighVaf = Doubles.greaterThan(depth.alleleReadCount(), 0.75 * depth.totalReadCount());
        if (isHighVaf) {
            return Doubles.lessThan(homPoisson(depth), 0.005) ? GermlineGenotypeStatus.HOM_ALT : GermlineGenotypeStatus.HET;
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
    private final String tumorSample;
    @NotNull
    private final Consumer<VariantContext> consumer;

    public GermlineGenotypeEnrichment(@NotNull final String germlineSample, @NotNull final String tumorSample,
            @NotNull final Consumer<VariantContext> consumer) {
        this.germlineSample = germlineSample;
        this.tumorSample = tumorSample;
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        Genotype germlineGenotype = context.getGenotype(germlineSample);
        Genotype tumorGenotype = context.getGenotype(tumorSample);
        AllelicDepth germlineDepth = AllelicDepth.fromGenotype(germlineGenotype);
        AllelicDepth tumorDepth = AllelicDepth.fromGenotype(tumorGenotype);

        GermlineGenotypeStatus germlineStatus = status(germlineDepth);
        GermlineGenotypeStatus tumorStatus = status(tumorDepth);
        GermlineGenotypeStatus status = combined(germlineStatus, tumorStatus);

        Allele refAllele = context.getReference();
        Allele altAllele = context.getAlternateAllele(0);
        List<Allele> genotypeAlleles = Lists.newArrayList();
        if (status.equals(GermlineGenotypeStatus.HOM_ALT)) {
            genotypeAlleles.add(altAllele);
        } else {
            genotypeAlleles.add(refAllele);
        }
        genotypeAlleles.add(altAllele);

        final Genotype updatedGenotype = new GenotypeBuilder(germlineGenotype).alleles(genotypeAlleles).make();
        final List<Genotype> updatedGenotypes = Lists.newArrayList();
        context.getGenotypes().stream().filter(x -> !x.getSampleName().equals(germlineSample)).forEach(updatedGenotypes::add);
        updatedGenotypes.add(updatedGenotype);

        VariantContextBuilder builder = new VariantContextBuilder(context).genotypes(updatedGenotypes);
        if (status.equals(GermlineGenotypeStatus.LOW_VAF)) {
            builder.filter(LOW_VAF_FILTER);
        }
        consumer.accept(builder.make());
    }

    @Override
    public void flush() {

    }

    static GermlineGenotypeStatus combined(GermlineGenotypeStatus germline, GermlineGenotypeStatus tumor) {
        if (germline == GermlineGenotypeStatus.HOM_ALT) {
            return tumor == GermlineGenotypeStatus.HOM_ALT ? GermlineGenotypeStatus.HOM_ALT : GermlineGenotypeStatus.HET;
        }

        return germline;
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
