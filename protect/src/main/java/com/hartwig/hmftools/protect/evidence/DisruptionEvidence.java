package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceType;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class DisruptionEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionableGenes;

    public DisruptionEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableGenes = actionableGenes.stream()
                .filter(x -> x.event() == GeneLevelEvent.ANY_MUTATION || x.event() == GeneLevelEvent.INACTIVATION
                        || x.event() == GeneLevelEvent.DELETION)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<ReportableHomozygousDisruption> reportables) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ReportableHomozygousDisruption reportable : reportables) {
            result.addAll(evidence(reportable));
        }
        return result;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull ReportableHomozygousDisruption reportable) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableGene actionable : actionableGenes) {
            if (actionable.gene().equals(reportable.gene())) {
                ProtectEvidence evidence = personalizedEvidenceFactory.somaticReportableEvidence(actionable)
                        .gene(reportable.gene())
                        .event("homozygous disruption")
                        .evidenceType(fromGeneLevelEvent(actionable.event()))
                        .build();
                result.add(evidence);
            }
        }

        return result;
    }

    @NotNull
    private static ProtectEvidenceType fromGeneLevelEvent(@NotNull GeneLevelEvent event) {
        switch (event) {
            case ANY_MUTATION:
                return ProtectEvidenceType.ANY_MUTATION;
            case DELETION:
                return ProtectEvidenceType.DELETION;
            case INACTIVATION:
                return ProtectEvidenceType.INACTIVATION;
            default:
                throw new IllegalStateException("Gene level event found in disruption evidence that should not exist: " + event);
        }
    }
}
