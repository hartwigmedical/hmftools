package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;

import org.jetbrains.annotations.NotNull;

public class CopyNumberEvidence {

    @NotNull
    private final List<ActionableGene> actionableGenes;

    public CopyNumberEvidence(@NotNull final List<ActionableGene> actionableGenes) {
        this.actionableGenes = actionableGenes;
    }

    @NotNull
    public List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull List<ReportableGainLoss> reportables) {
        List<ProtectEvidenceItem> result = Lists.newArrayList();
        for (ReportableGainLoss reportable : reportables) {
            result.addAll(evidence(doids, reportable));
        }
        return result;
    }

    @NotNull
    private List<ProtectEvidenceItem> evidence(@NotNull Set<String> doids, @NotNull ReportableGainLoss reportable) {
        List<ProtectEvidenceItem> result = Lists.newArrayList();
        for (ActionableGene actionable : actionableGenes) {
            if (actionable.gene().equals(reportable.gene()) && isTypeMatch(actionable, reportable)) {
                ProtectEvidenceItem evidence = ProtectEvidenceItems.builder(doids, actionable)
                        .germline(false)
                        .genomicEvent(reportable.genomicEvent())
                        .reported(true)
                        .build();
                result.add(evidence);
            }
        }
        return ProtectEvidenceItems.reportHighest(result);
    }

    private static boolean isTypeMatch(@NotNull ActionableGene actionable, @NotNull ReportableGainLoss reportable) {
        switch (actionable.event()) {
            case AMPLIFICATION:
                return reportable.interpretation() == CopyNumberInterpretation.GAIN;
            case INACTIVATION:
            case DELETION:
                return reportable.interpretation() == CopyNumberInterpretation.FULL_LOSS
                        || reportable.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS;
            default:
                return false;
        }
    }
}
