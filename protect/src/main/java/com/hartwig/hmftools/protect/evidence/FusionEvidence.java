package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class FusionEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionablePromiscuous;
    @NotNull
    private final List<ActionableFusion> actionableFusions;

    public FusionEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes, @NotNull final List<ActionableFusion> actionableFusions) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionablePromiscuous =
                actionableGenes.stream().filter(x -> x.event().equals(GeneLevelEvent.FUSION)).collect(Collectors.toList());
        this.actionableFusions = actionableFusions;
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<LinxFusion> fusions) {
        return fusions.stream().flatMap(x -> evidence(x).stream()).collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull LinxFusion reportable) {
        List<ProtectEvidence> geneEvidence = actionablePromiscuous.stream()
                .filter(x -> match(x, reportable))
                .map(x -> evidence(reportable, x))
                .collect(Collectors.toList());

        List<ProtectEvidence> fusionEvidence =
                actionableFusions.stream().filter(x -> match(x, reportable)).map(x -> evidence(reportable, x)).collect(Collectors.toList());

        List<ProtectEvidence> result = Lists.newArrayList();
        result.addAll(geneEvidence);
        result.addAll(fusionEvidence);

        return result;
    }

    @NotNull
    private ProtectEvidence evidence(@NotNull LinxFusion reportable, @NotNull ActionableEvent actionable) {
        return personalizedEvidenceFactory.somaticReportableEvidence(actionable).genomicEvent(reportable.genomicEvent()).build();
    }

    private static boolean match(@NotNull ActionableGene actionable, @NotNull LinxFusion reportable) {
        return actionable.gene().equals(reportable.geneStart()) || actionable.gene().equals(reportable.geneEnd());
    }

    private static boolean match(@NotNull ActionableFusion actionable, @NotNull LinxFusion reportable) {
        if (!actionable.geneDown().equals(reportable.geneEnd())) {
            return false;
        }

        if (!actionable.geneUp().equals(reportable.geneStart())) {
            return false;
        }

        Integer actionableMinExonDown = actionable.minExonDown();
        if (actionableMinExonDown != null && reportable.fusedExonDown() < actionableMinExonDown) {
            return false;
        }

        Integer actionableMaxExonDown = actionable.maxExonDown();
        if (actionableMaxExonDown != null && reportable.fusedExonDown() > actionableMaxExonDown) {
            return false;
        }

        Integer actionableMinExonUp = actionable.minExonUp();
        if (actionableMinExonUp != null && reportable.fusedExonUp() < actionableMinExonUp) {
            return false;
        }

        Integer actionableMaxExonUp = actionable.maxExonUp();
        if (actionableMaxExonUp != null && reportable.fusedExonUp() > actionableMaxExonUp) {
            return false;
        }

        return true;
    }
}
