package com.hartwig.hmftools.purple.somatic;

import static com.google.common.collect.Lists.newArrayList;

import java.io.File;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.clonality.ModifiableWeightedPloidy;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.clonality.PeakModelFactory;
import com.hartwig.hmftools.common.variant.enrich.SomaticPurityEnrichment;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SomaticPeakStream {

    private final SomaticFitConfig somaticFitConfig;
    private final CommonConfig commonConfig;
    private final String inputVCF;
    private final boolean enabled;

    private int indelCount;
    private int snpCount;

    public SomaticPeakStream(final ConfigSupplier configSupplier) {
        this.somaticFitConfig = configSupplier.somaticConfig();
        this.commonConfig = configSupplier.commonConfig();
        this.enabled = somaticFitConfig.file().isPresent();
        this.inputVCF = enabled ? somaticFitConfig.file().get().toString() : "";
    }

    public int indelCount() {
        return indelCount;
    }

    public int snpCount() {
        return snpCount;
    }

    public List<PeakModel> somaticPeakModel(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions) {

        if (enabled) {
            try (VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false)) {

                final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();
                final Consumer<VariantContext> consumer = context -> {
                    VariantContextDecorator decorator = new VariantContextDecorator(context);
                    if (Doubles.lessThan(decorator.variantCopyNumber(), somaticFitConfig.clonalityMaxPloidy()) && decorator.isPass()
                            && HumanChromosome.contains(decorator.chromosome()) && HumanChromosome.fromString(decorator.chromosome())
                            .isAutosome()) {
                        AllelicDepth depth = decorator.allelicDepth(commonConfig.tumorSample());
                        weightedPloidies.add(ModifiableWeightedPloidy.create()
                                .from(depth)
                                .setPloidy(decorator.variantCopyNumber())
                                .setWeight(1));
                    }

                    if (decorator.isPass()) {
                        if (decorator.type() == VariantType.INDEL) {
                            indelCount++;
                        } else {
                            snpCount++;
                        }
                    }
                };

                final SomaticPurityEnrichment somaticPurityEnrichment = new SomaticPurityEnrichment(commonConfig.version(),
                        commonConfig.tumorSample(),
                        purityAdjuster,
                        copyNumbers,
                        fittedRegions,
                        consumer);

                for (VariantContext context : vcfReader) {
                    somaticPurityEnrichment.accept(context);
                }

                return new PeakModelFactory(somaticFitConfig.clonalityMaxPloidy(), somaticFitConfig.clonalityBinWidth()).model(
                        weightedPloidies);
            }
        }

        return Lists.newArrayList();

    }
}
