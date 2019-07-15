package com.hartwig.hmftools.patientreporter.genepanel;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneModel {

    @NotNull
    public abstract Set<String> oncoDriverGenes();

    @NotNull
    public abstract Set<String> tsgDriverGenes();

    @NotNull
    public abstract Set<String> significantlyAmplifiedGenes();

    @NotNull
    public abstract Set<String> significantlyDeletedGenes();

    @Value.Derived
    public Set<String> somaticVariantGenes() {
        Set<String> somaticVariantGenePanel = Sets.newHashSet();
        somaticVariantGenePanel.addAll(oncoDriverGenes());
        somaticVariantGenePanel.addAll(tsgDriverGenes());
        return somaticVariantGenePanel;
    }

    @Value.Derived
    public boolean isDeletionReportable(@NotNull String gene) {
        return significantlyDeletedGenes().contains(gene) || tsgDriverGenes().contains(gene);
    }

    @Value.Derived
    public boolean isAmplificationReportable(@NotNull String gene) {
        return significantlyAmplifiedGenes().contains(gene) || oncoDriverGenes().contains(gene);
    }

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
