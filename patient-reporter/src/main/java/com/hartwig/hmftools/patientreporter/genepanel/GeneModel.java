package com.hartwig.hmftools.patientreporter.genepanel;

import java.util.Map;
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
    public abstract Set<String> somaticVariantDriverGenes();

    @NotNull
    public abstract Set<String> significantlyAmplifiedGenes();

    @NotNull
    public abstract Set<String> significantlyDeletedGenes();

    @NotNull
    public abstract Set<String> drupActionableGenes();

    @NotNull
    public abstract Map<String, DriverCategory> geneDriverCategoryMap();

    @Value.Derived
    public Set<String> somaticVariantGenes() {
        Set<String> somaticVariantGenePanel = Sets.newHashSet();
        somaticVariantGenePanel.addAll(somaticVariantDriverGenes());
        somaticVariantGenePanel.addAll(drupActionableGenes());
        return somaticVariantGenePanel;
    }

    @Value.Derived
    public boolean isDeletionReportable(@NotNull String gene) {
        return significantlyDeletedGenes().contains(gene) || geneDriverCategoryMap().get(gene) == DriverCategory.TSG;
    }

    @Value.Derived
    public boolean isAmplificationReportable(@NotNull String gene) {
        return significantlyAmplifiedGenes().contains(gene) || geneDriverCategoryMap().get(gene) == DriverCategory.ONCO;
    }
}
