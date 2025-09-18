package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.cobalt.exclusions.SuppliedExcludedRegions;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.gc.GCProfileFactory;

public class WindowStatuses implements CobaltCalculation.Filter
{
    private final ListMultimap<Chromosome, WindowStatus> mStatusesByChromosome = ArrayListMultimap.create();

    public WindowStatuses(CobaltConfig config) throws IOException
    {
        final ListMultimap<Chromosome, GCProfile> gcProfileData = GCProfileFactory.loadGCContent(WINDOW_SIZE, config.GcProfilePath);
        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(config.mExcludedRegions);
        ListMultimap<Chromosome, GCProfile> toExclude = excludedRegions.findIntersections(gcProfileData);
        gcProfileData.forEach( (chromosome, gcProfile) -> {
            boolean excluded = toExclude.containsEntry(chromosome, gcProfile);
            mStatusesByChromosome.put(chromosome, new WindowStatus(gcProfile, excluded));
        });
    }

    @Override
    public boolean exclude(final Chromosome chromosome, final ReadDepth readDepth)
    {
        List<WindowStatus> statusesForChromosome= mStatusesByChromosome.get(chromosome);
        WindowStatus status = statusesForChromosome.get(indexFor(readDepth));
        return status.maskedOut();
    }

    private int indexFor(ReadDepth depth)
    {
        return depth.StartPosition / WINDOW_SIZE;
    }
}
