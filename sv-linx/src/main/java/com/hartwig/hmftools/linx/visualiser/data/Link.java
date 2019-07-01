package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.types.ResolvedType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Link
{
    public abstract String sampleId();

    public abstract int clusterId();

    public abstract int chainId();

    public abstract int svId();

    public abstract StructuralVariantType type();

    public abstract ResolvedType resolvedType();

    public abstract String startChromosome();

    public abstract long startPosition();

    public abstract int startOrientation();

    public abstract String startInfo();

    public abstract String endChromosome();

    public abstract long endPosition();

    public abstract int endOrientation();

    public abstract String endInfo();

    public abstract int traverseCount();

    public boolean connectorsOnly()
    {
        return isSimpleSV() || isLineElement();
    }

    public boolean isSimpleSV()
    {
        return resolvedType().isSimple();
    }

    public boolean isLineElement()
    {
        return resolvedType() == ResolvedType.LINE;
    }

    public boolean isValidStart()
    {
        return HumanChromosome.contains(startChromosome());
    }

    public boolean isValidEnd()
    {
        return HumanChromosome.contains(endChromosome());
    }

}
