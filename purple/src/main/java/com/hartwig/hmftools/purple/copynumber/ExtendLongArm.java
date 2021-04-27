package com.hartwig.hmftools.purple.copynumber;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

final class ExtendLongArm
{
    private static final Set<HumanChromosome> FORCE_CHROMOSOMES = EnumSet.of(HumanChromosome._21);

    private static final Set<HumanChromosome> ELIGIBLE_CHROMOSOMES = EnumSet.of(HumanChromosome._13,
            HumanChromosome._14,
            HumanChromosome._15,
            HumanChromosome._21,
            HumanChromosome._22,
            HumanChromosome._Y);

    @NotNull
    static List<CombinedRegion> extendLongArm(@NotNull final List<CombinedRegion> regions)
    {
        final int centromereIndex = findCentromere(regions);
        if(centromereIndex > 0)
        {
            final double copyNumber = regions.get(centromereIndex).tumorCopyNumber();
            extendLeft(copyNumber, centromereIndex - 1, regions);
        }

        return regions;
    }

    private static int findCentromere(@NotNull final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            final String contig = regions.get(i).chromosome();
            if(!HumanChromosome.contains(contig) || !ELIGIBLE_CHROMOSOMES.contains(HumanChromosome.fromString(contig)))
            {
                return -1;
            }

            if(regions.get(i).region().support() == SegmentSupport.CENTROMERE)
            {
                return i;
            }
        }

        return -1;
    }

    private static boolean isProcessed(CombinedRegion region)
    {
        if(FORCE_CHROMOSOMES.contains(HumanChromosome.fromString(region.chromosome())))
        {
            return false;
        }

        return region.isProcessed();
    }

    private static void extendLeft(double copyNumber, int targetIndex, @NotNull final List<CombinedRegion> regions)
    {
        if(targetIndex < 0)
        {
            return;
        }

        final CombinedRegion target = regions.get(targetIndex);
        if(isProcessed(target))
        {
            return;
        }

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
