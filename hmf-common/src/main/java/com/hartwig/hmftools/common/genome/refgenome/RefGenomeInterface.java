package com.hartwig.hmftools.common.genome.refgenome;

import java.util.List;

public interface RefGenomeInterface
{
    String getBaseString(final String chromosome, int posStart, int posEnd);

    String getBaseString(final String chromosome, final List<int[]> baseRanges)
;}
