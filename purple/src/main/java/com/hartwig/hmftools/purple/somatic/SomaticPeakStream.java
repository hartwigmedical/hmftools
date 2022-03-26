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
import com.hartwig.hmftools.purple.SomaticVariantCache;
import com.hartwig.hmftools.purple.fitting.ModifiableWeightedPloidy;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.purple.fitting.PeakModelFactory;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SomaticPeakStream
{
    private final SomaticFitConfig mSomaticFitConfig;
    private final PurpleConfig mConfig;

    public SomaticPeakStream(final PurpleConfig config)
    {
        mSomaticFitConfig = config.SomaticFitting;
        mConfig = config;
    }

    public List<PeakModel> somaticPeakModel(final SomaticVariantCache somaticVariants)
    {
        if(!somaticVariants.hasData())
            return Lists.newArrayList();

        final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();

        for(SomaticData variant : somaticVariants.variants())
        {
            if(!variant.isPass())
                continue;

            if(variant.variantCopyNumber() >= mSomaticFitConfig.clonalityMaxPloidy())
                continue;

            if(!HumanChromosome.contains(variant.chromosome()) || !HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            AllelicDepth depth = variant.tumorAlleleDepth();

            if(depth != null)
            {
                weightedPloidies.add(ModifiableWeightedPloidy.create()
                        .from(depth)
                        .setPloidy(variant.variantCopyNumber())
                        .setWeight(1));
            }
        }

        return new PeakModelFactory(
                mSomaticFitConfig.clonalityMaxPloidy(), mSomaticFitConfig.clonalityBinWidth()).model(weightedPloidies);
    }
}
