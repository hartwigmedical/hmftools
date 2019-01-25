package com.hartwig.hmftools.patientreporter.genepanel;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneModel {

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> somaticVariantDriverGenePanel();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> significantlyAmplifiedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> significantlyDeletedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> drupActionableGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> disruptionGeneWhiteList();

    @NotNull
    public abstract Map<String, DriverCategory> geneDriverCategoryMap();

    @Value.Derived
    public Set<String> somaticVariantGenePanel() {
        Set<String> somaticVariantGenePanel = Sets.newHashSet();
        somaticVariantGenePanel.addAll(somaticVariantDriverGenePanel().keySet());
        somaticVariantGenePanel.addAll(drupActionableGenes().keySet());
        return somaticVariantGenePanel;
    }

    @Value.Derived
    public boolean isDeletionReportable(@NotNull String gene) {
        return significantlyDeletedGenes().keySet().contains(gene) || geneDriverCategoryMap().get(gene) == DriverCategory.TSG;
    }

    @Value.Derived
    public boolean isAmplificationReportable(@NotNull String gene) {
        return significantlyAmplifiedGenes().keySet().contains(gene) || geneDriverCategoryMap().get(gene) == DriverCategory.ONCO;
    }

    @Value.Derived
    @NotNull
    public Set<String> disruptionGenePanel() {
        Set<String> disruptionGenePanel = Sets.newHashSet();
        for (String driverGene : somaticVariantDriverGenePanel().keySet()) {
            if (geneDriverCategoryMap().get(driverGene) != DriverCategory.ONCO) {
                disruptionGenePanel.add(driverGene);
            }
        }

        for (String drupActionableGene : drupActionableGenes().keySet()) {
            if (geneDriverCategoryMap().get(drupActionableGene) != DriverCategory.ONCO) {
                disruptionGenePanel.add(drupActionableGene);
            }
        }

        disruptionGenePanel.addAll(disruptionGeneWhiteList().keySet());

        return disruptionGenePanel;
    }
}
