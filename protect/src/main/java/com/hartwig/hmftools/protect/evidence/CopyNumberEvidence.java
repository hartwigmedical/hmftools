package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.EventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class CopyNumberEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionableGenes;

    public CopyNumberEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableGenes = actionableGenes.stream()
                .filter(x -> x.event() == GeneLevelEvent.INACTIVATION || x.event() == GeneLevelEvent.AMPLIFICATION
                        || x.event() == GeneLevelEvent.OVEREXPRESSION || x.event() == GeneLevelEvent.DELETION
                        || x.event() == GeneLevelEvent.UNDEREXPRESSION || x.event() == GeneLevelEvent.ANY_MUTATION)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<GainLoss> reportableGainsLosses, @NotNull List<GainLoss> allGainsLosses) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (GainLoss reportableGainLoss : reportableGainsLosses) {
            result.addAll(evidence(reportableGainLoss, true));
        }

        for (GainLoss gainLoss : allGainsLosses) {
            if (!reportableGainsLosses.contains(gainLoss)) {
                result.addAll(evidence(gainLoss, false));
            }
        }

        return result;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull GainLoss gainLoss, boolean report) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableGene actionable : actionableGenes) {
            if (actionable.gene().equals(gainLoss.gene()) && isTypeMatch(actionable, gainLoss)) {
                ProtectEvidence evidence = personalizedEvidenceFactory.somaticEvidence(actionable)
                        .reported(report)
                        .gene(gainLoss.gene())
                        .transcript(gainLoss.transcript())
                        .isCanonical(gainLoss.isCanonical())
                        .event(EventGenerator.copyNumberEvent(gainLoss))
                        .eventIsHighDriver(EvidenceDriverLikelihood.interpretCopyNumber())
                        .build();
                result.add(evidence);
            }
        }

        return result;
    }

    private static boolean isTypeMatch(@NotNull ActionableGene actionable, @NotNull GainLoss reportable) {
        switch (actionable.event()) {
            case AMPLIFICATION:
            case OVEREXPRESSION:
                return reportable.interpretation() == CopyNumberInterpretation.FULL_GAIN
                        || reportable.interpretation() == CopyNumberInterpretation.PARTIAL_GAIN;
            case INACTIVATION:
            case DELETION:
            case UNDEREXPRESSION:
            case ANY_MUTATION:
                return reportable.interpretation() == CopyNumberInterpretation.FULL_LOSS
                        || reportable.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS;
            default:
                throw new IllegalStateException(
                        "Actionable event found in copy number evidence that should not exist: " + actionable.event());
        }
    }
}
