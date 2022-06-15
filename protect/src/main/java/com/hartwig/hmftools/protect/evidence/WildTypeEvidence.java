package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public class WildTypeEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionableGenes;
    @NotNull
    private final Map<String, DriverGene> driverGenes;

    public WildTypeEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes, @NotNull final Map<String, DriverGene> driverGenes) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionableGenes = actionableGenes.stream().filter(x -> x.event() == GeneLevelEvent.WILD_TYPE).collect(Collectors.toList());
        this.driverGenes = driverGenes;
    }

    public List<ProtectEvidence> evidence(@NotNull List<ReportableVariant> reportableGermlineVariant,
            @NotNull List<ReportableVariant> reportableSomaticVariant, @NotNull List<GainLoss> reportableSomaticGainsLosses,
            @NotNull List<LinxFusion> reportableFusions, @NotNull List<ReportableHomozygousDisruption> homozygousDisruptions) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (ActionableGene wildType : actionableGenes) {

            DriverGene driverGene = driverGenes.get(wildType.gene());

            if (driverGene != null) {
                boolean wildTypeSomaticVariant = false;
                boolean wildTypeGermlineVariant = false;
                boolean wildTypeSomaticGainLoss = false;
                boolean wildTypeFusions = false;
                boolean wildTypeHomozygousDisruption = false;

                List<Boolean> somaticVariantBoolean = Lists.newArrayList();
                List<Boolean> somaticGermlineBoolean = Lists.newArrayList();
                List<Boolean> CNVBoolean = Lists.newArrayList();
                List<Boolean> FusionBooleanEnd = Lists.newArrayList();
                List<Boolean> FusionBooleanStart = Lists.newArrayList();
                List<Boolean> HomozygousDisruptionBoolean = Lists.newArrayList();

                for (ReportableVariant somaticVariant : reportableSomaticVariant) {
                    isWildType(driverGene.gene(), somaticVariant.gene(), somaticVariantBoolean);
                    wildTypeSomaticVariant = somaticVariantBoolean.contains(true);
                }

                for (ReportableVariant germlineVariant : reportableGermlineVariant) {
                    isWildType(driverGene.gene(), germlineVariant.gene(), somaticGermlineBoolean);
                    wildTypeGermlineVariant = somaticGermlineBoolean.contains(true);
                }

                for (GainLoss gainLoss : reportableSomaticGainsLosses) {
                    isWildType(driverGene.gene(), gainLoss.gene(), CNVBoolean);
                    wildTypeSomaticGainLoss = CNVBoolean.contains(true);
                }

                for (LinxFusion fusions : reportableFusions) {
                    isWildType(driverGene.gene(), fusions.geneStart(), FusionBooleanStart);
                    isWildType(driverGene.gene(), fusions.geneEnd(), FusionBooleanEnd);
                    wildTypeFusions =
                            (!FusionBooleanStart.contains(true) && !FusionBooleanEnd.contains(true)) || (FusionBooleanStart.contains(true)
                                    && FusionBooleanEnd.contains(true));
                }

                for (ReportableHomozygousDisruption homozygousDisruption : homozygousDisruptions) {

                    isWildType(driverGene.gene(), homozygousDisruption.gene(), HomozygousDisruptionBoolean);
                    wildTypeHomozygousDisruption = HomozygousDisruptionBoolean.contains(true);
                }

                if (wildTypeSomaticVariant && wildTypeGermlineVariant && wildTypeSomaticGainLoss && wildTypeFusions
                        && wildTypeHomozygousDisruption) {
                    evidences.add(evidence(wildType));
                }
            }
        }
        return evidences;
    }

    public void isWildType(@NotNull String driverGene, @NotNull String reportableGene, @NotNull List<Boolean> booleanList) {
        if (driverGene.equals(reportableGene)) {
            booleanList.add(false);
        } else {
            booleanList.add(true);
        }
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