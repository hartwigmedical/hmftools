package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public interface Variant
{
    CategoryType categoryType();

    String description();

    String gene();

    String sequence();

    default List<String> refSequences() { return Lists.newArrayList(); }

    double copyNumber();

    double vaf();

    int tumorFragments();

    default boolean hasPhaseVariants() { return false; }

    boolean reported();

    void generateSequences(final RefGenomeInterface refGenome, final PvConfig config);

}
