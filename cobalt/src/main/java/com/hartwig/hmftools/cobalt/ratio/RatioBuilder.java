package com.hartwig.hmftools.cobalt.ratio;

import com.google.common.collect.ArrayListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

public interface RatioBuilder
{
    ArrayListMultimap<Chromosome, ReadRatio> ratios();
}
