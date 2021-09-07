package com.hartwig.hmftools.common.purple.cnchromosome;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CnPerChromosomeArmData {

    @NotNull
    public abstract HumanChromosome chromosome();

    @NotNull
    public abstract ChromosomeArm chromosomeArm();

    public abstract double copyNumber();

}
