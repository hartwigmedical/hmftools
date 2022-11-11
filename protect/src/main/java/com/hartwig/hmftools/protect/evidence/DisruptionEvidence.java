package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.serve.datamodel.gene.ActionableGene;
import com.hartwig.serve.datamodel.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class DisruptionEvidence {

    static final String HOMOZYGOUS_DISRUPTION_EVENT = "homozygous disruption";

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionableGenes;

    public DisruptionEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableGenes = actionableGenes.stream()
                .filter(x -> x.event() == GeneLevelEvent.ANY_MUTATION || x.event() == GeneLevelEvent.INACTIVATION
                        || x.event() == GeneLevelEvent.DELETION || x.event() == GeneLevelEvent.UNDEREXPRESSION)
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<HomozygousDisruption> reportables) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (HomozygousDisruption reportable : reportables) {
            result.addAll(evidence(reportable));
        }
        return result;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull HomozygousDisruption reportable) {
        List<ProtectEvidence> result = Lists.newArrayList();
        for (ActionableGene actionable : actionableGenes) {
            if (actionable.gene().equals(reportable.gene())) {
                ProtectEvidence evidence = personalizedEvidenceFactory.somaticReportableEvidence(actionable)
                        .gene(reportable.gene())
                        .transcript(reportable.transcript())
                        .isCanonical(reportable.isCanonical())
                        .event(HOMOZYGOUS_DISRUPTION_EVENT)
                        .eventIsHighDriver(true)
                        .build();
                result.add(evidence);
            }
        }

        return result;
    }
}
