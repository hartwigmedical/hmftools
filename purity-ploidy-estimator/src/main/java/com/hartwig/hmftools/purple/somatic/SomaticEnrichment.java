package com.hartwig.hmftools.purple.somatic;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SomaticEnrichment {

    private static final String PURPLE_CN_INFO = "PURPLE_CN";
    private static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_INFO = "PURPLE_MAP";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";
    private static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    private static final String PURPLE_AF_INFO = "PURPLE_AF";
    private static final String PURPLE_AF_DESC = "Purity adjusted allelic frequency of variant";
    private static final String PURPLE_PLOIDY_INFO = "PURPLE_PLOIDY";
    private static final String PURPLE_PLOIDY_DESC = "Purity adjusted ploidy of variant";

    @NotNull
    private final String tumorSample;
    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;
    @NotNull
    private final GenomeRegionSelector<FittedRegion> fittedRegionSelector;

    SomaticEnrichment(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions, final String tumorSample) {
        this(tumorSample, purityAdjuster, Multimaps.fromRegions(copyNumbers), Multimaps.fromRegions(fittedRegions));
    }

    private SomaticEnrichment(@NotNull final String tumorSample, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<Chromosome, PurpleCopyNumber> copyNumbers,
            @NotNull final Multimap<Chromosome, FittedRegion> fittedRegions) {
        this.tumorSample = tumorSample;
        this.purityAdjuster = purityAdjuster;
        this.copyNumberSelector = GenomeRegionSelectorFactory.createImproved(copyNumbers);
        this.fittedRegionSelector = GenomeRegionSelectorFactory.createImproved(fittedRegions);
    }

    @NotNull
    @VisibleForTesting
    static VCFHeader generateOutputHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_PLOIDY_INFO, 1, VCFHeaderLineType.Float, PURPLE_PLOIDY_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_PLOIDY_INFO,
                1,
                VCFHeaderLineType.Float,
                PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        return outputVCFHeader;
    }

    @NotNull
    public VariantContext enrich(@NotNull final VariantContext variant) {

        final Genotype genotype = variant.getGenotype(tumorSample);
        if (genotype != null && genotype.hasAD() && HumanChromosome.contains(variant.getContig())) {
            final VariantContextBuilder builder = new VariantContextBuilder(variant);

            final GenomePosition position = GenomePositions.create(variant.getContig(), variant.getStart());
            fittedRegionSelector.select(position).ifPresent(x -> enrich(x, builder));
            copyNumberSelector.select(position).ifPresent(x -> enrich(x, SomaticVariantFactory.allelicDepth(genotype), builder));

            return builder.make();
        }
        return variant;
    }

    @NotNull
    private VariantContextBuilder enrich(@NotNull final FittedRegion fittedRegion, @NotNull final VariantContextBuilder builder) {
        return builder.attribute(PURPLE_GERMLINE_INFO, fittedRegion.status().toString());
    }

    @NotNull
    private VariantContextBuilder enrich(@NotNull final PurpleCopyNumber purpleCopyNumber, @NotNull final AllelicDepth depth,
            @NotNull final VariantContextBuilder builder) {

        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double vaf = purityAdjuster.purityAdjustedVAF(purpleCopyNumber.chromosome(), Math.max(0.001, copyNumber), depth.alleleFrequency());

        return builder.attribute(PURPLE_CN_INFO, copyNumber)
                .attribute(PURPLE_PLOIDY_INFO, vaf * copyNumber)
                .attribute(PURPLE_AF_INFO, vaf)
                .attribute(PURPLE_MINOR_ALLELE_PLOIDY_INFO, purpleCopyNumber.minorAllelePloidy());

    }

}
