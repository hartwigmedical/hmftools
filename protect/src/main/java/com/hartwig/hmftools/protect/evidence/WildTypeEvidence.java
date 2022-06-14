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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class WildTypeEvidence {
    private static final Logger LOGGER = LogManager.getLogger(WildTypeEvidence.class);

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
            LOGGER.info(wildType);

            DriverGene driverGene = driverGenes.get(wildType.gene());

            if (driverGene != null) {
                LOGGER.info(driverGene);

                boolean wildTypeSomaticVariant = false;
                boolean wildTypeGermlineVariant = false;
                boolean wildTypeSomaticGainLoss = false;
                boolean wildTypeFusions = false;
                boolean wildTypeHomozygousDisruption = false;

                for (ReportableVariant somaticVariant: reportableSomaticVariant) {
                    wildTypeSomaticVariant = determineWildTypes(driverGene.gene(), somaticVariant.gene());
                }

                for (ReportableVariant germlineVariant: reportableGermlineVariant) {
                    wildTypeGermlineVariant = determineWildTypes(driverGene.gene(), germlineVariant.gene());
                }

                for (GainLoss gainLoss: reportableSomaticGainsLosses){
                    wildTypeSomaticGainLoss = determineWildTypes(driverGene.gene(), gainLoss.gene());
                }

                for (LinxFusion fusions: reportableFusions){
                    boolean wildTypeFusionsStart = determineWildTypes(driverGene.gene(), fusions.geneStart());
                    boolean wildTypeFusionsEnd = determineWildTypes(driverGene.gene(), fusions.geneEnd());
                    wildTypeFusions = !wildTypeFusionsStart && !wildTypeFusionsEnd;
                }

                for (ReportableHomozygousDisruption homozygousDisruption: homozygousDisruptions){
                    wildTypeHomozygousDisruption = determineWildTypes(driverGene.gene(), homozygousDisruption.gene());
                }

                LOGGER.info("wildTypeSomaticVariant: " +wildTypeSomaticVariant);
                LOGGER.info("wildTypeGermlineVariant: " + wildTypeGermlineVariant);
                LOGGER.info("wildTypeSomaticGainLoss: " + wildTypeSomaticGainLoss);
                LOGGER.info("wildTypeFusions: " + wildTypeFusions);
                LOGGER.info("wildTypeHomozygousDisruption: " + wildTypeHomozygousDisruption);
                if (wildTypeSomaticVariant && wildTypeGermlineVariant && wildTypeSomaticGainLoss && wildTypeFusions
                        && wildTypeHomozygousDisruption) {
                    evidences.add(evidence(wildType));
                }
            }
        }
        LOGGER.info(evidences);
        return evidences;
    }

    public boolean determineWildTypes(@NotNull String driverGene, @NotNull String reportableGene) {
        if (driverGene.equals(reportableGene)) {
            return false;
        } else {
            return true;
        }
    }

    @NotNull
    private ProtectEvidence evidence( @NotNull ActionableGene actionable) {
        return personalizedEvidenceFactory.somaticEvidence(actionable)
                .reported(false)
                .gene(actionable.gene())
                .event(actionable.gene() + " wild type")
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretWildType())
                .build();
    }
}