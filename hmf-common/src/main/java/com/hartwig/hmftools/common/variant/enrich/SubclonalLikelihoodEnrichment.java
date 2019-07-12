package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.clonality.PeakModelFactory;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihood;
import com.hartwig.hmftools.common.variant.clonality.SubclonalLikelihoodFactory;
import com.hartwig.hmftools.common.variant.filter.PassingFilter;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SubclonalLikelihoodEnrichment implements VariantContextEnrichment {

    private static final Logger LOGGER = LogManager.getLogger(SubclonalLikelihoodEnrichment.class);

    private static final double MAX_PLOIDY = 10;
    private static final double BIN_WIDTH = 0.05;

    public static final String SUBCLONAL_LIKELIHOOD_FLAG = "SUBCL";
    private static final String SUBCLONAL_LIKELIHOOD_FLAG_DESCIPTION = "Subclonal likelihood";

    private final String sample;
    private final Consumer<VariantContext> consumer;

    private final List<VariantContext> buffer = Lists.newArrayList();
    private final List<ModifiableWeightedPloidy> ploidies = Lists.newArrayList();
    private final PeakModelFactory peakModelFactory;
    private final List<PeakModel> peakModel = Lists.newArrayList();
    private final VariantContextFilter filter;

    public SubclonalLikelihoodEnrichment(final String sample, @NotNull final Consumer<VariantContext> consumer) {
        this.sample = sample;
        this.consumer = consumer;
        this.peakModelFactory = new PeakModelFactory(MAX_PLOIDY, BIN_WIDTH);
        this.filter = new PassingFilter();
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        buffer.add(context);

        final Genotype genotype = context.getGenotype(sample);

        if (genotype != null && genotype.hasAD() && context.hasAttribute(PurityEnrichment.PURPLE_PLOIDY_INFO) && filter.test(context)
                && HumanChromosome.contains(context.getContig()) && HumanChromosome.fromString(context.getContig()).isAutosome()) {

            double ploidy = context.getAttributeAsDouble(PurityEnrichment.PURPLE_PLOIDY_INFO, MAX_PLOIDY);
            if (Doubles.lessThan(ploidy, MAX_PLOIDY)) {
                final AllelicDepth depth = AllelicDepth.fromGenotype(genotype);
                final ModifiableWeightedPloidy weightedPloidy = ModifiableWeightedPloidy.create()
                        .setPloidy(ploidy)
                        .setAlleleReadCount(depth.alleleReadCount())
                        .setTotalReadCount(depth.totalReadCount())
                        .setWeight(1)
                        .setDistribution(new BinomialDistribution(depth.totalReadCount(), depth.alleleFrequency()));

                ploidies.add(weightedPloidy);
            }
        }
    }

    @NotNull
    public List<PeakModel> peakModel() {
        return peakModel;
    }

    @Override
    public void flush() {
        LOGGER.info("Enriching somatics with subclonal likelihood");
        peakModel.addAll(peakModelFactory.model(ploidies));
        final List<SubclonalLikelihood> likelihoods = SubclonalLikelihoodFactory.subclonalLikelihood(BIN_WIDTH, peakModel);

        double[] subclonalLikelihoods = new double[(int) Math.round(MAX_PLOIDY / BIN_WIDTH)];
        for (SubclonalLikelihood likelihood : likelihoods) {
            if (!Doubles.isZero(likelihood.likelihood())) {
                int key = (int) Math.round(likelihood.bucket() / BIN_WIDTH);
                subclonalLikelihoods[key] = likelihood.likelihood();
            }
        }

        for (final VariantContext context : buffer) {

            final double ploidy = context.getAttributeAsDouble(PurityEnrichment.PURPLE_PLOIDY_INFO, MAX_PLOIDY);
            if (Doubles.lessThan(ploidy, MAX_PLOIDY)) {
                int key = (int) Math.round(ploidy / BIN_WIDTH);
                double variantLikelihood = Math.round(subclonalLikelihoods[key] * 1000d) / 1000d;
                if (!Doubles.isZero(variantLikelihood)) {
                    context.getCommonInfo().putAttribute(SUBCLONAL_LIKELIHOOD_FLAG, variantLikelihood);
                }
                consumer.accept(context);
            }
        }

        buffer.clear();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFInfoHeaderLine(SUBCLONAL_LIKELIHOOD_FLAG,
                1,
                VCFHeaderLineType.Float,
                SUBCLONAL_LIKELIHOOD_FLAG_DESCIPTION));

        return template;
    }
}
