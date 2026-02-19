package com.hartwig.hmftools.amber.blacklist;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

public record AmberBlacklistPoint(String chromosome, int position, int count, double meanVaf) implements GenomePosition
{
}
