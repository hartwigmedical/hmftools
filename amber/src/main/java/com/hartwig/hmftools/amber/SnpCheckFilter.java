package com.hartwig.hmftools.amber;

import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

public class SnpCheckFilter implements Predicate<PositionEvidence>
{
    private final Set<GenomePosition> mSnpLoci;

    public SnpCheckFilter(final Multimap<Chromosome, AmberSite> snpLoci)
    {
        mSnpLoci = snpLoci.values()
                .stream()
                .filter(AmberSite::snpCheck)
                .map(x -> GenomePositions.create(x.chromosome(), x.position()))
                .collect(Collectors.toSet());
    }

    @Override
    public boolean test(final PositionEvidence baseDepth)
    {
        return mSnpLoci.contains(GenomePositions.create(baseDepth));
    }
}
