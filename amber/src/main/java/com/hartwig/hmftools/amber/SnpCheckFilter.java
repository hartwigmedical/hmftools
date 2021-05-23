package com.hartwig.hmftools.amber;

import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class SnpCheckFilter implements Predicate<BaseDepth>
{
    private final Set<GenomePosition> mSnpLoci;

    public SnpCheckFilter(@NotNull final Multimap<Chromosome, AmberSite> snpLoci)
    {
        mSnpLoci = snpLoci.values()
                .stream()
                .filter(AmberSite::snpCheck)
                .map(x -> GenomePositions.create(x.chromosome(), x.position()))
                .collect(Collectors.toSet());
    }

    @Override
    public boolean test(final BaseDepth baseDepth)
    {
        return mSnpLoci.contains(GenomePositions.create(baseDepth));
    }
}
