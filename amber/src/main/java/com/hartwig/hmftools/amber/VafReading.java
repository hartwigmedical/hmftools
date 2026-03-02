package com.hartwig.hmftools.amber;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public record VafReading(String chromosome, int position, int readDepth, int refSupport, int altSupport, double minorAlleleFrequency)
        implements GenomePosition
{
    public VafReading(String chromosome, int position, int readDepth, int refSupport, int altSupport)
    {
        this(chromosome, position, readDepth, refSupport, altSupport, 0.5);
    }
}
