package com.hartwig.hmftools.ctdna;

import java.util.List;
import java.util.Map;

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

    double gc();

    int tumorFragments();

    default boolean hasPhaseVariants() { return false; }

    boolean reported();

    void generateSequences(final RefGenomeInterface refGenome, final PvConfig config);

    boolean checkAndRegisterLocation(final Map<String,List<Integer>> registeredLocations);

    default int sequenceCount() { return 1 + refSequences().size(); }
}
