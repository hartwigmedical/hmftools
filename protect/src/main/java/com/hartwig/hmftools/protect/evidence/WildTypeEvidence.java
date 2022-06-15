package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.wildtype.WildTypeGene;
import com.hartwig.hmftools.common.wildtype.WildTypeFactory;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class WildTypeEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionableGenes;
    @NotNull
    private final List<DriverGene> driverGenes;

    public WildTypeEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes, @NotNull final List<DriverGene> driverGenes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableGenes = actionableGenes.stream().filter(x -> x.event() == GeneLevelEvent.WILD_TYPE).collect(Collectors.toList());
        this.driverGenes = driverGenes;
    }

    public List<ProtectEvidence> evidence(@NotNull List<ReportableVariant> reportableGermlineVariant,
            @NotNull List<ReportableVariant> reportableSomaticVariant, @NotNull List<GainLoss> reportableSomaticGainsLosses,
            @NotNull List<LinxFusion> reportableFusions, @NotNull List<ReportableHomozygousDisruption> homozygousDisruptions,
            @NotNull List<ReportableGeneDisruption> geneDisruptions) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        List<WildTypeGene> wildTypeGenes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                geneDisruptions,
                driverGenes);

        for (ActionableGene actionable : actionableGenes) {
            for (WildTypeGene wildType : wildTypeGenes) {
                if (wildType.gene().equals(actionable.gene())) {
                    evidences.add(evidence(actionable));
                }
            }
        }
        return evidences;
    }

    @NotNull
    private ProtectEvidence evidence(@NotNull ActionableGene actionable) {
        return personalizedEvidenceFactory.somaticEvidence(actionable)
                .reported(false)
                .gene(actionable.gene())
                .event(actionable.gene() + " wild type")
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretWildType())
                .build();
    }
}