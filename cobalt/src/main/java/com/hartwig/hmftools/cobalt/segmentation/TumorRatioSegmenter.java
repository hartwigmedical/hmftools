package com.hartwig.hmftools.cobalt.segmentation;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

class TumorRatioSegmenter extends CobaltRatioSegmenter
{
    TumorRatioSegmenter(final ListMultimap<Chromosome, CobaltRatio> ratios, final ChrArmLocator chrArmLocator, final double gamma)
    {
        super(ratios, chrArmLocator, gamma);
    }

    @Override
    public double value(final CobaltRatio ratio)
    {
        return ratio.tumorGCRatio();
    }
}
