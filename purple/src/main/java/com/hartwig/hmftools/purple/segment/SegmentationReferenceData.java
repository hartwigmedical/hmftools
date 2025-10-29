package com.hartwig.hmftools.purple.segment;

import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.purple.ReferenceData;

public class SegmentationReferenceData
{
    private final ReferenceData mReferenceData;

    public SegmentationReferenceData(final ReferenceData referenceData)
    {
        this.mReferenceData = referenceData;
    }

    public Map<Chromosome, GenomePosition> chromosomeLengths()
    {
        return mReferenceData.ChromosomeLengths;
    }

    public Map<Chromosome, GenomePosition> centromeres()
    {
        return mReferenceData.Centromeres;
    }
}
