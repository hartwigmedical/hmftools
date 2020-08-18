package com.hartwig.hmftools.common.drivercatalog.panel;

import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DriverGenePanel {

    @NotNull
    public abstract Map<String, DndsDriverGeneLikelihood> tsgLikelihood();

    @NotNull
    public abstract Map<String, DndsDriverImpactLikelihood> oncoLikelihood();

    @NotNull
    public abstract Set<String> amplificationTargets();

    @NotNull
    public abstract Set<String> deletionTargets();

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

}
