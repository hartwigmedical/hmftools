package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.clonality.PeakModelFactory;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihood;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment implements VariantContextEnrichment {

    public static final String SUBCLONAL_LIKELIHOOD_FLAG = "SUBCL";
    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION = "Non-zero subclonal likelihood";

    private final String tumorSample;
    private final double clonalityMaxPloidy;
    private final double clonalityBinWidth;
    private final Consumer<VariantContext> consumer;
    private final List<VariantContextDecorator> buffer = Lists.newArrayList();
    private final List<PeakModel> peakModel = Lists.newArrayList();

    SubclonalLikelihoodEnrichment(String tumorSample, double maxPloidy, double binWidth, @NotNull final Consumer<VariantContext> consumer) {
        this.clonalityMaxPloidy = maxPloidy;
        this.clonalityBinWidth = binWidth;
        this.consumer = consumer;
        this.tumorSample = tumorSample;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        buffer.add(new VariantContextDecorator(context));
    }

    @Override
    public void flush() {
        peakModel.clear();
        peakModel.addAll(createPeakModel());
        SubclonalLikelihood likelihoodFactory = new SubclonalLikelihood(clonalityBinWidth, peakModel);

        for (VariantContextDecorator variant : buffer) {
            double subclonalLikelihood = Math.round(likelihoodFactory.subclonalLikelihood(variant.variantCopyNumber()) * 1000d) / 1000d;

            if (!Doubles.isZero(subclonalLikelihood)) {
                variant.context().getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, subclonalLikelihood);
            }
            consumer.accept(variant.context());
        }

        buffer.clear();
    }

    @NotNull
    public List<PeakModel> somaticPeakModel() {
        return peakModel;
    }

    @NotNull
    private List<PeakModel> createPeakModel() {
        final List<ModifiableWeightedPloidy> weightedPloidies = Lists.newArrayList();

        for (VariantContextDecorator context : buffer) {
            if (Doubles.lessThan(context.variantCopyNumber(), clonalityMaxPloidy) && context.isPass()
                    && HumanChromosome.contains(context.chromosome()) && HumanChromosome.fromString(context.chromosome()).isAutosome()) {
                AllelicDepth depth = context.allelicDepth(tumorSample);
                weightedPloidies.add(ModifiableWeightedPloidy.create().from(depth).setPloidy(context.variantCopyNumber()).setWeight(1));
            }
        }
        return new PeakModelFactory(clonalityMaxPloidy, clonalityBinWidth).model(weightedPloidies);
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,
                1,
                VCFHeaderLineType.Float,
                SUBCLONAL_LIKELIHOOD_FLAG_DESCRIPTION));

        return template;
    }
}
