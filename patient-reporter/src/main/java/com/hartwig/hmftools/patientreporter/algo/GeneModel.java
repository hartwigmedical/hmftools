package com.hartwig.hmftools.patientreporter.algo;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.region.GenomeRegion;
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
    public abstract Map<String, DriverCategory> somaticVariantDriverCategoryMap();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> significantlyAmplifiedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> significantlyDeletedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> drupActionableGenes();

    @Value.Derived
    @Nullable
    public DriverCategory geneDriverCategory(@NotNull String gene) {
        return somaticVariantDriverCategoryMap().get(gene);
    }

    @Value.Derived
    public Set<String> somaticVariantGenePanel() {
        Set<String> somaticVariantGenePanel = Sets.newHashSet();
        somaticVariantGenePanel.addAll(somaticVariantDriverGenePanel().keySet());
        somaticVariantGenePanel.addAll(drupActionableGenes().keySet());
        return somaticVariantGenePanel;
    }

    @Value.Derived
    @NotNull
    public Set<String> cnvGenePanel() {
        Set<String> cnvGenes = Sets.newHashSet();
        cnvGenes.addAll(somaticVariantDriverGenePanel().keySet());
        cnvGenes.addAll(significantlyAmplifiedGenes().keySet());
        cnvGenes.addAll(significantlyDeletedGenes().keySet());
        cnvGenes.addAll(drupActionableGenes().keySet());
        return cnvGenes;
    }

    @Value.Derived
    public boolean isDeletionReportable(@NotNull String gene) {
        // TODO (KODU): Filter out onco genes from drup actionable genes.
        return significantlyDeletedGenes().keySet().contains(gene) || somaticVariantDriverCategoryMap().get(gene) != DriverCategory.ONCO
                || drupActionableGenes().keySet().contains(gene);
    }

    @Value.Derived
    public boolean isAmplificationReportable(@NotNull String gene) {
        // TODO (KODU): Filter out tsg genes from drup actionable genes.
        return significantlyAmplifiedGenes().keySet().contains(gene) || somaticVariantDriverCategoryMap().get(gene) != DriverCategory.TSG
                || drupActionableGenes().keySet().contains(gene);
    }

    @Value.Derived
    @NotNull
    public Set<String> disruptionGeneIDPanel() {
        // KODU: Structural variant analyser requires a set of ensembl IDs rather than a set of gene names.
        Set<String> disruptionGeneIDPanel = Sets.newHashSet();
        for (Map.Entry<String, HmfTranscriptRegion> driverGene : somaticVariantDriverGenePanel().entrySet()) {
            if (somaticVariantDriverCategoryMap().get(driverGene.getKey()) == DriverCategory.TSG) {
                disruptionGeneIDPanel.add(driverGene.getValue().geneID());
            }
        }

        for (Map.Entry<String, HmfTranscriptRegion> drupActionableGene : drupActionableGenes().entrySet()) {
            if (somaticVariantDriverCategoryMap().get(drupActionableGene.getKey()) != DriverCategory.ONCO) {
                disruptionGeneIDPanel.add(drupActionableGene.getValue().geneID());
            }
        }

        return disruptionGeneIDPanel;
    }

    @Value.Derived
    public long somaticVariantsNumberOfBases() {
        return somaticVariantDriverGenePanel().values().stream().mapToLong(GenomeRegion::bases).sum();
    }

    @Value.Derived
    public int somaticVariantNumberOfGenes() {
        return somaticVariantDriverGenePanel().size();
    }
}
