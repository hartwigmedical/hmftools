package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DriverGenePanel
{
    @NotNull
    public abstract List<DriverGene> driverGenes();

    @NotNull
    public abstract Map<String, DndsDriverGeneLikelihood> tsgLikelihood();

    @NotNull
    public abstract Map<String, DndsDriverGeneLikelihood> oncoLikelihood();

    @NotNull
    public Set<DriverGene> amplificationTargets()
    {
        return targets(DriverGene::reportAmplification);
    }

    @NotNull
    public Set<DriverGene> deletionTargets()
    {
        return targets(DriverGene::reportDeletion);
    }

    @NotNull
    private Set<DriverGene> targets(@NotNull Predicate<DriverGene> filter)
    {
        return driverGenes().stream().filter(filter).collect(Collectors.toSet());
    }
}
