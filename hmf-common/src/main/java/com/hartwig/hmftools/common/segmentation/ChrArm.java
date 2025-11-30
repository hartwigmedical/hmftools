package com.hartwig.hmftools.common.segmentation;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public record ChrArm(HumanChromosome chromosome, Arm arm) implements Comparable<ChrArm>
{
    @Override
    public int compareTo(@NotNull final ChrArm o)
    {
        int cmp = chromosome.compareTo(o.chromosome);
        if(cmp == 0)
        {
            return arm.compareTo(o.arm);
        }
        return cmp;
    }
}
