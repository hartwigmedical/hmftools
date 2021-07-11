package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
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
    public List<ProtectEvidence> evidence(@NotNull List<LinxFusion> reportableFusions, @NotNull List<LinxFusion> unreportedFusions) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (LinxFusion reportable : reportableFusions) {
            evidences.addAll(evidence(reportable));
        }

        for (LinxFusion unreported : unreportedFusions) {
            evidences.addAll(evidence(unreported));
        }
        return evidences;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull LinxFusion fusion) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (ActionableGene promiscuous : actionablePromiscuous) {
            if (match(fusion, promiscuous)) {
                evidences.add(evidence(fusion, promiscuous));
            }
        }

        for (ActionableFusion actionableFusion : actionableFusions) {
            if (match(fusion, actionableFusion)) {
                evidences.add(evidence(fusion, actionableFusion));
            }
        }
        return evidences;
    }

    @NotNull
    private ProtectEvidence evidence(@NotNull LinxFusion fusion, @NotNull ActionableEvent actionable) {
        return personalizedEvidenceFactory.somaticEvidence(actionable)
                .reported(fusion.reported())
                .genomicEvent(fusion.genomicEvent())
                .build();
    }

    private static boolean match(@NotNull LinxFusion fusion, @NotNull ActionableGene actionable) {
        return actionable.gene().equals(fusion.geneStart()) || actionable.gene().equals(fusion.geneEnd());
    }

    private static boolean match(@NotNull LinxFusion fusion, @NotNull ActionableFusion actionable) {
        if (!actionable.geneDown().equals(fusion.geneEnd())) {
            return false;
        }

        if (!actionable.geneUp().equals(fusion.geneStart())) {
            return false;
        }

        Integer actionableMinExonDown = actionable.minExonDown();
        if (actionableMinExonDown != null && fusion.fusedExonDown() < actionableMinExonDown) {
            return false;
        }

        Integer actionableMaxExonDown = actionable.maxExonDown();
        if (actionableMaxExonDown != null && fusion.fusedExonDown() > actionableMaxExonDown) {
            return false;
        }

        Integer actionableMinExonUp = actionable.minExonUp();
        if (actionableMinExonUp != null && fusion.fusedExonUp() < actionableMinExonUp) {
            return false;
        }

        Integer actionableMaxExonUp = actionable.maxExonUp();
        if (actionableMaxExonUp != null && fusion.fusedExonUp() > actionableMaxExonUp) {
            return false;
        }

        return true;
    }
}
