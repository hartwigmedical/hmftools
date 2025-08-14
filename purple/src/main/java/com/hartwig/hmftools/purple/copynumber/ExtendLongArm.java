package com.hartwig.hmftools.purple.copynumber;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.SegmentSupport;

public final class ExtendLongArm
{
    private static final Set<HumanChromosome> FORCE_CHROMOSOMES = Sets.newHashSet();
    private static final Set<HumanChromosome> ELIGIBLE_CHROMOSOMES = Sets.newHashSet();

    static
    {
        Arrays.stream(HumanChromosome.values()).filter(x -> x.hasShortArm()).forEach(x -> FORCE_CHROMOSOMES.add(x));
        ELIGIBLE_CHROMOSOMES.addAll(FORCE_CHROMOSOMES);
        ELIGIBLE_CHROMOSOMES.add(HumanChromosome._Y);
    }

    public static List<CombinedRegion> extendLongArm(final List<CombinedRegion> regions)
    {
        final int centromereIndex = findCentromere(regions);
        if(centromereIndex > 0)
        {
            CombinedRegion centromere = regions.get(centromereIndex);

            if(centromere.copyNumberMethod() != CopyNumberMethod.UNKNOWN)
            {
                double copyNumber = regions.get(centromereIndex).tumorCopyNumber();
                extendLeft(copyNumber, centromereIndex - 1, regions);
            }
        }

        return regions;
    }

    private static int findCentromere(final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            final String contig = regions.get(i).chromosome();

            if(!HumanChromosome.contains(contig) || !ELIGIBLE_CHROMOSOMES.contains(HumanChromosome.fromString(contig)))
                return -1;

            if(regions.get(i).region().support() == SegmentSupport.CENTROMERE)
                return i;
        }

        return -1;
    }

    private static boolean isProcessed(CombinedRegion region)
    {
        if(FORCE_CHROMOSOMES.contains(HumanChromosome.fromString(region.chromosome())))
            return false;

        return region.isProcessed();
    }

    private static void extendLeft(double copyNumber, int targetIndex, final List<CombinedRegion> regions)
    {
        if(targetIndex < 0)
            return;

        final CombinedRegion target = regions.get(targetIndex);

        if(isProcessed(target))
            return;

        target.setTumorCopyNumber(CopyNumberMethod.LONG_ARM, copyNumber);
        if(target.region().support() != SegmentSupport.NONE)
        {
            extendLeft(copyNumber, targetIndex - 1, regions);
            return;
        }

        targetIndex--;
        while(targetIndex >= 0)
        {
            final CombinedRegion neighbour = regions.get(targetIndex);
            if(isProcessed(neighbour))
            {
                return;
            }

            target.extend(neighbour.region());
            regions.remove(targetIndex);
            targetIndex--;

            if(target.region().support() != SegmentSupport.NONE)
            {
                extendLeft(copyNumber, targetIndex, regions);
                return;
            }
        }
    }
}
