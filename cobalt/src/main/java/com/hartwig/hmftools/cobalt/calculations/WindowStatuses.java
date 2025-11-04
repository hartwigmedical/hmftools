package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.WINDOW_SIZE;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.cobalt.diploid.DiploidStatus;
import com.hartwig.hmftools.cobalt.exclusions.SuppliedExcludedRegions;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

class WindowStatuses
{
    private final ListMultimap<Chromosome, WindowStatus> mStatusesByChromosome = ArrayListMultimap.create();
    private final ListMultimap<Chromosome, GCProfile> GcProfileData;

    WindowStatuses(ListMultimap<Chromosome, GCProfile> gcProfileData,
            List<ChrBaseRegion> exclusions,
            ListMultimap<Chromosome, DiploidStatus> diploidRegions)
    {
        GcProfileData = gcProfileData;
        boolean checkDiploid = !diploidRegions.isEmpty();
        SuppliedExcludedRegions excludedRegions = new SuppliedExcludedRegions(exclusions);
        ListMultimap<Chromosome, GCProfile> toExclude = excludedRegions.findIntersections(gcProfileData);
        gcProfileData.keySet().forEach(chromosome ->
        {
            List<DiploidStatus> diploidStatuses = diploidRegions.get(chromosome);
            List<GCProfile> gcProfiles = gcProfileData.get(chromosome);
            gcProfiles.forEach(gcProfile ->
            {
                boolean excluded = toExclude.containsEntry(chromosome, gcProfile);
                boolean nonDiploid = false;
                if(checkDiploid)
                {
                    int index = indexFor(gcProfile.start());
                    if(index < diploidStatuses.size())
                    {
                        nonDiploid = !diploidStatuses.get(index).isDiploid;
                    }
                    else
                    {
                        nonDiploid = true;
                    }
                }
                mStatusesByChromosome.put(chromosome, new WindowStatus(gcProfile, excluded, nonDiploid));
            });
        });
    }

    Double referenceGcValueForWindow(final Chromosome chromosome, int position)
    {
        return GcProfileData.get(chromosome).get(indexFor(position)).gcContent();
    }

    boolean exclude(final Chromosome chromosome, final DepthReading readDepth)
    {
        List<WindowStatus> statusesForChromosome = mStatusesByChromosome.get(chromosome);
        WindowStatus status = statusesForChromosome.get(indexFor(readDepth.StartPosition));
        return status.maskedOut();
    }

    private int indexFor(int position)
    {
        return position / WINDOW_SIZE;
    }
}
