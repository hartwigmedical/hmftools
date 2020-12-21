package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class FusionEvidence {

    @NotNull
    private final List<ActionableGene> actionableGenes;
    @NotNull
    private final List<ActionableFusion> actionableFusions;

    public FusionEvidence(@NotNull final List<ActionableGene> actionableGenes, @NotNull final List<ActionableFusion> actionableFusions) {
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

        return ProtectEvidenceItems.reportHighest(result);
    }

    private static boolean match(@NotNull ActionableGene actionable, @NotNull LinxFusion victim) {
        return actionable.gene().equals(victim.geneStart()) || actionable.gene().equals(victim.geneEnd());
    }

    private static boolean match(@NotNull ActionableFusion actionable, @NotNull LinxFusion victim) {
        if (!actionable.geneDown().equals(victim.geneStart())) {
            return false;
        }

        if (!actionable.geneUp().equals(victim.geneEnd())) {
            return false;
        }

        Integer actionableMinExonDown = actionable.minExonDown();
        if (actionableMinExonDown != null && victim.fusedExonDown() < actionableMinExonDown) {
            return false;
        }

        Integer actionableMaxExonDown = actionable.maxExonDown();
        if (actionableMaxExonDown != null && victim.fusedExonDown() > actionableMaxExonDown) {
            return false;
        }

        Integer actionableMinExonUp = actionable.minExonUp();
        if (actionableMinExonUp != null && victim.fusedExonUp() < actionableMinExonUp) {
            return false;
        }

        Integer actionableMaxExonUp = actionable.maxExonUp();
        if (actionableMaxExonUp != null && victim.fusedExonUp() > actionableMaxExonUp) {
            return false;
        }

        return true;
    }

    @NotNull
    private static ProtectEvidenceItem evidence(@NotNull Set<String> doid, @NotNull LinxFusion reportable,
            @NotNull ActionableEvent actionable) {
        return ProtectEvidenceItems.builder(doid, actionable).genomicEvent(reportable.genomicEvent()).reported(true).build();
    }
}
