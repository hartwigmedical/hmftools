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
import com.hartwig.hmftools.purple.fitting.ModifiableWeightedPloidy;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.purple.fitting.PeakModelFactory;
import com.hartwig.hmftools.common.variant.enrich.SomaticPurityEnrichment;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SomaticPeakStream
{
    private final SomaticFitConfig mSomaticFitConfig;
    private final CommonConfig mCommonConfig;
    private final String mSomaticVcf;
    private final boolean mEnabled;

    private int mIndelCount;
    private int mSnpCount;

    public SomaticPeakStream(final ConfigSupplier configSupplier)
    {
        mSomaticFitConfig = configSupplier.somaticConfig();
        mCommonConfig = configSupplier.commonConfig();
        mEnabled = mSomaticFitConfig.file().isPresent();
        mSomaticVcf = mEnabled ? mSomaticFitConfig.file().get().toString() : "";
    }

    public int indelCount()
    {
        return mIndelCount;
    }

    public int snpCount()
    {
        return mSnpCount;
    }

    public List<PeakModel> somaticPeakModel(
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions)
    {
        if(!mEnabled)
            return Lists.newArrayList();

        try (VCFFileReader vcfReader = new VCFFileReader(new File(mSomaticVcf), false))
        {
            // gather up passing variants less then max ploidy
            final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();

            final Consumer<VariantContext> consumer = context ->
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(Doubles.lessThan(decorator.variantCopyNumber(), mSomaticFitConfig.clonalityMaxPloidy())
                && decorator.isPass()
                && HumanChromosome.contains(decorator.chromosome()) && HumanChromosome.fromString(decorator.chromosome()).isAutosome())
                {
                    AllelicDepth depth = decorator.allelicDepth(mCommonConfig.tumorSample());
                    weightedPloidies.add(ModifiableWeightedPloidy.create()
                            .from(depth)
                            .setPloidy(decorator.variantCopyNumber())
                            .setWeight(1));
                }

                if(decorator.isPass())
                {
                    if(decorator.type() == VariantType.INDEL)
                    {
                        mIndelCount++;
                    }
                    else
                    {
                        mSnpCount++;
                    }
                }
            };

            final SomaticPurityEnrichment somaticPurityEnrichment = new SomaticPurityEnrichment(
                    mCommonConfig.version(), mCommonConfig.tumorSample(),
                    purityAdjuster, copyNumbers, fittedRegions, consumer);

            for(VariantContext context : vcfReader)
            {
                somaticPurityEnrichment.accept(context);
            }

            return new PeakModelFactory(
                    mSomaticFitConfig.clonalityMaxPloidy(), mSomaticFitConfig.clonalityBinWidth()).model(weightedPloidies);
        }
    }
}
