package com.hartwig.hmftools.protect.common;

import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DriverGeneView {

    @NotNull
    public abstract Set<String> oncoDriverGenes();

    @NotNull
    public abstract Set<String> tsgDriverGenes();

    @Value.Derived
    @Nullable
    public DriverCategory category(@NotNull String gene) {
        if (oncoDriverGenes().contains(gene)) {
            return DriverCategory.ONCO;
        } else if (tsgDriverGenes().contains(gene)) {
            return DriverCategory.TSG;
        }
        return null;
    }
}