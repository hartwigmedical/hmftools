package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionEvidence {

    private final List<ActionableGene> actionableGenes;
    private final List<ActionableFusion> actionableFusions;

    public FusionEvidence(final List<ActionableGene> actionableGenes, final List<ActionableFusion> actionableFusions) {
        this.actionableGenes = actionableGenes.stream().filter(x -> x.event().equals(GeneLevelEvent.FUSION)).collect(Collectors.toList());
        this.actionableFusions = actionableFusions;
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doid, @NotNull List<LinxFusion> fusions) {
        return fusions.stream().flatMap(x -> evidence(doid, x).stream()).collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doid, @NotNull LinxFusion reportable) {
        List<ProtectEvidenceItem> geneEvidence = actionableGenes.stream()
                .filter(x -> match(x, reportable))
                .map(x -> evidence(doid, reportable, x))
                .collect(Collectors.toList());

        List<ProtectEvidenceItem> fusionEvidence = actionableFusions.stream()
                .filter(x -> match(x, reportable))
                .map(x -> evidence(doid, reportable, x))
                .collect(Collectors.toList());

        Set<ProtectEvidenceItem> result = Sets.newHashSet();
        result.addAll(geneEvidence);
        result.addAll(fusionEvidence);

        return ProtectEvidenceItems.doNotReportInsignificantEvidence(result);
    }

    private static boolean match(ActionableGene actionable, LinxFusion victim) {
        return actionable.gene().equals(victim.geneStart()) || actionable.gene().equals(victim.geneEnd());
    }

    private static boolean match(ActionableFusion actionable, LinxFusion victim) {
        if (!actionable.geneDown().equals(victim.geneStart())) {
            return false;
        }

        if (!actionable.geneUp().equals(victim.geneEnd())) {
            return false;
        }

        @Nullable
        Integer actionableExonDown = actionable.exonDown();
        if (actionableExonDown != null && actionableExonDown != victim.fusedExonDown()) {
            return false;
        }

        @Nullable
        Integer actionableExonUp = actionable.exonUp();
        return actionableExonUp == null || actionableExonUp == victim.fusedExonUp();
    }

    @NotNull
    private static ProtectEvidenceItem evidence(@NotNull Set<String> doid, @NotNull LinxFusion reportable,
            @NotNull ActionableEvent actionable) {
        return ProtectEvidenceItems.builder(doid, actionable).genomicEvent(reportable.genomicEvent()).reported(true).build();
    }

}
