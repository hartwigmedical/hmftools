package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DriverGenePanel {

    @NotNull
    public abstract List<DriverGene> driverGenes();

    @NotNull
    public abstract Map<String, DndsDriverGeneLikelihood> tsgLikelihood();

    @NotNull
    public abstract Map<String, DndsDriverGeneLikelihood> oncoLikelihood();

    @NotNull
    public Set<String> amplificationTargets() {
        return targets(DriverGene::reportAmplification);
    }

    @NotNull
    public Set<String> deletionTargets() {
        return targets(DriverGene::reportDeletion);
    }

    @NotNull
    public Set<String> disruptionTargets() {
        return targets(DriverGene::reportDisruption);
    }

    @NotNull
    public abstract Map<String, String> deletionBandMap();

    @NotNull
    public Set<String> oncoGenes() {
        return oncoLikelihood().keySet();
    }

    @NotNull
    public Set<String> tsGenes() {
        return tsgLikelihood().keySet();
    }

    @Nullable
    public DriverCategory category(@NotNull final String gene) {
        if (oncoLikelihood().containsKey(gene)) {
            return DriverCategory.ONCO;
        } else if (tsgLikelihood().containsKey(gene)) {
            return DriverCategory.TSG;
        }
        return null;
    }

    private Set<String> targets(Predicate<DriverGene> filter) {
        return driverGenes().stream().filter(filter).map(DriverGene::gene).collect(Collectors.toSet());
    }

}
