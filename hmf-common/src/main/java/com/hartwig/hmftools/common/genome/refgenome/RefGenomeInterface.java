package com.hartwig.hmftools.common.genome.refgenome;

import java.util.List;
import java.util.Map;

public interface RefGenomeInterface
{
    String getBaseString(final String chromosome, int posStart, int posEnd);

    String getBaseString(final String chromosome, final List<int[]> baseRanges);

    int getChromosomeLength(final String chromosome);

    byte[] getBases(final String chromosome, int posStart, int posEnd);

    default String getBase(final String chromosome, int pos) { return getBaseString(chromosome, pos, pos); }

    Map<String,Integer> chromosomeLengths();

    default boolean oneBasedIndexing() { return true; }
}
