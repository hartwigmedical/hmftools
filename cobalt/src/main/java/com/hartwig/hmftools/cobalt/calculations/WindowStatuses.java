package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.exclusions.SuppliedExcludedRegions;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class WindowStatuses implements GenomeFilter
{
    private final ListMultimap<Chromosome, WindowStatus> mStatusesByChromosome = ArrayListMultimap.create();

    public WindowStatuses(ListMultimap<Chromosome, GCProfile> gcProfileData, List<ChrBaseRegion> exclusions)
    {
        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(exclusions);
        ListMultimap<Chromosome, GCProfile> toExclude = excludedRegions.findIntersections(gcProfileData);
        gcProfileData.forEach( (chromosome, gcProfile) -> {
            boolean excluded = toExclude.containsEntry(chromosome, gcProfile);
            mStatusesByChromosome.put(chromosome, new WindowStatus(gcProfile, excluded));
        });
    }

    @Override
    public boolean exclude(final Chromosome chromosome, final DepthReading readDepth)
    {
        List<WindowStatus> statusesForChromosome= mStatusesByChromosome.get(chromosome);
        WindowStatus status = statusesForChromosome.get(indexFor(readDepth));
        return status.maskedOut();
    }

    private int indexFor(DepthReading depth)
    {
        return depth.StartPosition / WINDOW_SIZE;
    }
}
