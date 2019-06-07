package com.hartwig.hmftools.amber;

import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

class SnpCheckFilter implements Predicate<BaseDepth> {

    private final Set<GenomePosition> snpLoci;

    SnpCheckFilter(@NotNull final SortedSetMultimap<String, GenomeRegion> snpLoci) {
        this.snpLoci = snpLoci.values().stream().map(x -> GenomePositions.create(x.chromosome(), x.start())).collect(Collectors.toSet());
    }

    @Override
    public boolean test(final BaseDepth baseDepth) {
        return snpLoci.contains(GenomePositions.create(baseDepth));
    }
}
