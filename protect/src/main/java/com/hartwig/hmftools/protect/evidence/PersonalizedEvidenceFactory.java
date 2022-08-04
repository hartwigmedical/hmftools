package com.hartwig.hmftools.protect.evidence;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.ImmutableKnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.KnowledgebaseSource;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.cancertype.CancerType;
import com.hartwig.hmftools.serve.cancertype.CancerTypeFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PersonalizedEvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(PersonalizedEvidenceFactory.class);

    @NotNull
    private final Set<String> patientTumorDoids;

    @NotNull
    private final DoidParents doidParentModel;

    public PersonalizedEvidenceFactory(@NotNull final Set<String> patientTumorDoids, @NotNull final DoidParents doidParentModel) {
        this.patientTumorDoids = patientTumorDoids;
        this.doidParentModel = doidParentModel;
    }

    @NotNull
    public ImmutableProtectEvidence.Builder somaticEvidence(@NotNull ActionableEvent event) {
        return evidenceBuilder(event).reported(false).germline(false);
    }

    @NotNull
    public ImmutableProtectEvidence.Builder somaticReportableEvidence(@NotNull ActionableEvent event) {
        return evidenceBuilder(event).reported(true).germline(false);
    }

    @NotNull
    public ImmutableProtectEvidence.Builder evidenceBuilder(@NotNull ActionableEvent actionable) {
        return ImmutableProtectEvidence.builder()
                .treatment(actionable.treatment().treament())
                .onLabel(isOnLabel(actionable.applicableCancerType(), actionable.blacklistCancerTypes(), actionable.treatment().treament()))
                .level(actionable.level())
                .direction(actionable.direction())
                .sources(Sets.newHashSet(resolveProtectSource(actionable)));
    }

    public boolean isOnLabel(@NotNull CancerType applicableCancerType, @NotNull Set<CancerType> blacklistCancerTypes,
            @NotNull String treatment) {
        return patientTumorDoids.contains(applicableCancerType.doid()) && !isBlacklisted(blacklistCancerTypes, treatment);
    }

    @VisibleForTesting
    boolean isBlacklisted(@NotNull Set<CancerType> blacklistCancerTypes, @NotNull String treatment) {
        Set<String> blacklistDoids = CancerTypeFactory.doidStrings(blacklistCancerTypes);
        Set<String> allDoids = Sets.newHashSet();

        if (!blacklistDoids.isEmpty()) {
            LOGGER.info(" Starting doid resolving for blacklisting evidence  '{}' for treatment '{}'", blacklistDoids, treatment);
        }

        for (String doid : blacklistDoids) {
            allDoids.add(doid);
            allDoids.addAll(doidParentModel.parents(doid));
        }

        for (String result : allDoids) {
            for (String doidPatient : patientTumorDoids) {
                if (doidPatient.equals(result)) {
                    return true;
                }
            }
        }
        return false;
    }

    @NotNull
    private static KnowledgebaseSource resolveProtectSource(@NotNull ActionableEvent actionable) {
        return ImmutableKnowledgebaseSource.builder()
                .name(actionable.source())
                .sourceEvent(actionable.sourceEvent())
                .sourceUrls(actionable.sourceUrls())
                .evidenceType(determineEvidenceType(actionable))
                .rangeRank(determineRangeRank(actionable))
                .evidenceUrls(actionable.evidenceUrls())
                .build();
    }

    @VisibleForTesting
    @Nullable
    static Integer determineRangeRank(@NotNull ActionableEvent actionable) {
        return actionable instanceof ActionableRange ? ((ActionableRange) actionable).rank() : null;
    }

    @VisibleForTesting
    @NotNull
    static EvidenceType determineEvidenceType(@NotNull ActionableEvent actionable) {
        if (actionable instanceof ActionableHotspot) {
            return EvidenceType.HOTSPOT_MUTATION;
        } else if (actionable instanceof ActionableRange) {
            return fromActionableRange((ActionableRange) actionable);
        } else if (actionable instanceof ActionableGene) {
            return fromActionableGene((ActionableGene) actionable);
        } else if (actionable instanceof ActionableFusion) {
            return EvidenceType.FUSION_PAIR;
        } else if (actionable instanceof ActionableCharacteristic) {
            return fromActionableCharacteristic((ActionableCharacteristic) actionable);
        } else if (actionable instanceof ActionableHLA) {
            return EvidenceType.HLA;
        } else {
            throw new IllegalStateException("Unexpected actionable event detected in variant evidence: " + actionable);
        }
    }

    @NotNull
    private static EvidenceType fromActionableRange(@NotNull ActionableRange range) {
        switch (range.rangeType()) {
            case EXON:
                return EvidenceType.EXON_MUTATION;
            case CODON:
                return EvidenceType.CODON_MUTATION;
            default: {
                throw new IllegalStateException("Unsupported range type: " + range.rangeType());
            }
        }
    }

    @NotNull
    private static EvidenceType fromActionableGene(@NotNull ActionableGene gene) {
        switch (gene.event()) {
            case AMPLIFICATION:
                return EvidenceType.AMPLIFICATION;
            case OVEREXPRESSION:
                return EvidenceType.OVER_EXPRESSION;
            case DELETION:
                return EvidenceType.DELETION;
            case UNDEREXPRESSION:
                return EvidenceType.UNDER_EXPRESSION;
            case ACTIVATION:
                return EvidenceType.ACTIVATION;
            case INACTIVATION:
                return EvidenceType.INACTIVATION;
            case ANY_MUTATION:
                return EvidenceType.ANY_MUTATION;
            case FUSION:
                return EvidenceType.PROMISCUOUS_FUSION;
            case WILD_TYPE:
                return EvidenceType.WILD_TYPE;
            default: {
                throw new IllegalStateException("Unsupported gene level event: " + gene.event());
            }
        }
    }

    @NotNull
    private static EvidenceType fromActionableCharacteristic(@NotNull ActionableCharacteristic characteristic) {
        switch (characteristic.name()) {
            case MICROSATELLITE_UNSTABLE:
            case MICROSATELLITE_STABLE:
            case HIGH_TUMOR_MUTATIONAL_LOAD:
            case LOW_TUMOR_MUTATIONAL_LOAD:
            case HIGH_TUMOR_MUTATIONAL_BURDEN:
            case LOW_TUMOR_MUTATIONAL_BURDEN:
            case HOMOLOGOUS_RECOMBINATION_DEFICIENT:
                return EvidenceType.SIGNATURE;
            case HPV_POSITIVE:
            case EBV_POSITIVE:
                return EvidenceType.VIRAL_PRESENCE;
            default: {
                throw new IllegalStateException("Unsupported tumor characteristic: " + characteristic.name());
            }
        }
    }
}