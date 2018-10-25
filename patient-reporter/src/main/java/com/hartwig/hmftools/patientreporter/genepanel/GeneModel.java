package com.hartwig.hmftools.patientreporter.genepanel;

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
    public abstract Map<String, HmfTranscriptRegion> significantlyAmplifiedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> significantlyDeletedGenes();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> drupActionableGenes();

    @NotNull
    public abstract Map<String, DriverCategory> geneDriverCategoryMap();

    @Value.Derived
    @Nullable
    public DriverCategory geneDriverCategory(@NotNull String gene) {
        return geneDriverCategoryMap().get(gene);
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
        return significantlyDeletedGenes().keySet().contains(gene) || geneDriverCategoryMap().get(gene) != DriverCategory.ONCO;
    }

    @Value.Derived
    public boolean isAmplificationReportable(@NotNull String gene) {
        return significantlyAmplifiedGenes().keySet().contains(gene) || geneDriverCategoryMap().get(gene) != DriverCategory.TSG;
    }

    @Value.Derived
    @NotNull
    public Set<String> disruptionGeneIDPanel() {
        // KODU: Structural variant analyser requires a set of ensembl IDs rather than a set of gene names.
        Set<String> disruptionGeneIDPanel = Sets.newHashSet();
        for (Map.Entry<String, HmfTranscriptRegion> driverGene : somaticVariantDriverGenePanel().entrySet()) {
            if (geneDriverCategoryMap().get(driverGene.getKey()) != DriverCategory.ONCO) {
                disruptionGeneIDPanel.add(driverGene.getValue().geneID());
            }
        }

        for (Map.Entry<String, HmfTranscriptRegion> drupActionableGene : drupActionableGenes().entrySet()) {
            if (geneDriverCategoryMap().get(drupActionableGene.getKey()) != DriverCategory.ONCO) {
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
