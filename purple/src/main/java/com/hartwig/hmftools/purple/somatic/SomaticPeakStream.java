package com.hartwig.hmftools.purple.somatic;

import static java.lang.String.format;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CLONALITY_BIN_WIDTH;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CLONALITY_MAX_PLOIDY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.purple.fitting.ModifiableWeightedPloidy;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.purple.fitting.PeakModelFactory;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.SomaticFitConfig;

public class SomaticPeakStream
{
    private final SomaticFitConfig mSomaticFitConfig;

    public SomaticPeakStream(final PurpleConfig config)
    {
        mSomaticFitConfig = config.SomaticFitting;
    }

    public List<PeakModel> somaticPeakModel(final SomaticVariantCache somaticVariants)
    {
        if(!somaticVariants.hasData())
            return Lists.newArrayList();

        final List<ModifiableWeightedPloidy> weightedPloidies = newArrayList();

        for(SomaticVariant variant : somaticVariants.variants())
        {
            if(!variant.isPass())
                continue;

            if(variant.copyNumber() >= CLONALITY_MAX_PLOIDY)
                continue;

            if(!HumanChromosome.contains(variant.chromosome()) || !HumanChromosome.fromString(variant.chromosome()).isAutosome())
                continue;

            AllelicDepth depth = variant.tumorAlleleDepth();

            if(depth != null)
            {
                weightedPloidies.add(ModifiableWeightedPloidy.create()
                        .from(depth)
                        .setPloidy(variant.copyNumber())
                        .setWeight(1));
            }
        }

        PPL_LOGGER.debug("somatic peak uses {} variants", weightedPloidies.size());

        PeakModelFactory modelFactory = new PeakModelFactory(CLONALITY_MAX_PLOIDY, CLONALITY_BIN_WIDTH);

        List<PeakModel> peakModels = modelFactory.model(weightedPloidies);

        for(PeakModel peakModel : peakModels)
        {
            PPL_LOGGER.trace(format("somatic peak(%.4f wgt=%.4f) bucket(%.4f wgt=%.4f) valid(%s) subclonal(%s)",
                    peakModel.peak(), peakModel.peakAvgWeight(), peakModel.bucket(), peakModel.bucketWeight(),
                    peakModel.isValid(), peakModel.isSubclonal()));
        }

        return peakModels;
    }
}
