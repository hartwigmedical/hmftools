package com.hartwig.hmftools.protect.evidence;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.protect.EventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FusionEvidence {

    @NotNull
    private final PersonalizedEvidenceFactory personalizedEvidenceFactory;
    @NotNull
    private final List<ActionableGene> actionablePromiscuous;
    @NotNull
    private final List<ActionableFusion> actionableFusions;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public FusionEvidence(@NotNull final PersonalizedEvidenceFactory personalizedEvidenceFactory,
            @NotNull final List<ActionableGene> actionableGenes, @NotNull final List<ActionableFusion> actionableFusions,
            @NotNull final KnownFusionCache fusionCache) {
        this.personalizedEvidenceFactory = personalizedEvidenceFactory;
        this.actionablePromiscuous = actionableGenes.stream()
                .filter(x -> x.event().equals(GeneLevelEvent.FUSION) || x.event().equals(GeneLevelEvent.ACTIVATION) || x.event()
                        .equals(GeneLevelEvent.ANY_MUTATION))
                .collect(Collectors.toList());
        this.actionableFusions = actionableFusions;
        this.knownFusionCache = fusionCache;
    }

    @NotNull
    public List<ProtectEvidence> evidence(@NotNull List<LinxFusion> reportableFusions, @NotNull List<LinxFusion> allFusions) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (LinxFusion reportable : reportableFusions) {
            evidences.addAll(evidence(reportable));
        }

        for (LinxFusion allFusion : allFusions) {
            if (!allFusion.reported()) {
                evidences.addAll(evidence(allFusion));
            }
        }

        return evidences;
    }

    @NotNull
    private List<ProtectEvidence> evidence(@NotNull LinxFusion fusion) {
        List<ProtectEvidence> evidences = Lists.newArrayList();
        for (ActionableGene promiscuous : actionablePromiscuous) {
            if (promiscuous.event().equals(GeneLevelEvent.FUSION) && match(fusion, promiscuous)) {
                evidences.add(evidence(fusion, promiscuous));
            }
            if ((promiscuous.event().equals(GeneLevelEvent.ACTIVATION) || promiscuous.event().equals(GeneLevelEvent.ANY_MUTATION))
                    && match(fusion, promiscuous)) {
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
                .gene(geneFromActionable(actionable))
                .event(EventGenerator.fusionEvent(fusion))
                .eventIsHighDriver(EvidenceDriverLikelihood.interpretFusion(fusion.likelihood()))
                .build();
    }

    private boolean match(@NotNull LinxFusion fusion, @NotNull ActionableGene actionable) {
        if (fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_3.name())) {
            if (knownFusionCache.hasPromiscuousThreeGene(fusion.geneStart())) {
                return actionable.gene().equals(fusion.geneStart());
            }
        } else if (fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString())) {
            if (knownFusionCache.hasPromiscuousFiveGene(fusion.geneEnd())) {
                return actionable.gene().equals(fusion.geneEnd());
            }
        } else if (fusion.reportedType().equals(KnownFusionType.IG_PROMISCUOUS.toString())) {
            if (knownFusionCache.hasPromiscuousIgFusion(fusion.geneStart())) {
                return actionable.gene().equals(fusion.geneStart());
            } else if (knownFusionCache.hasPromiscuousIgFusion(fusion.geneEnd())) {
                return actionable.gene().equals(fusion.geneEnd());
            }
        } else if (fusion.reportedType().equals(KnownFusionType.NONE.toString())) {
            if (knownFusionCache.hasPromiscuousThreeGene(fusion.geneStart())) {
                return actionable.gene().equals(fusion.geneStart());
            } else if (knownFusionCache.hasPromiscuousFiveGene(fusion.geneEnd())) {
                return actionable.gene().equals(fusion.geneEnd());
            } else {
                return actionable.gene().equals(fusion.geneEnd()) || actionable.gene().equals(fusion.geneStart());
            }
        }
        return false;
    }

    private static boolean match(@NotNull LinxFusion fusion, @NotNull ActionableFusion actionable) {
        if (fusion.reportedType().equals(KnownFusionType.KNOWN_PAIR.toString()) || fusion.reportedType()
                .equals(KnownFusionType.EXON_DEL_DUP.toString()) || fusion.reportedType()
                .equals(KnownFusionType.IG_KNOWN_PAIR.toString())) {
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
        return false;
    }

    @Nullable
    private static String geneFromActionable(@NotNull ActionableEvent actionable) {
        if (actionable instanceof ActionableGene) {
            return ((ActionableGene) actionable).gene();
        } else if (actionable instanceof ActionableFusion) {
            return null;
        } else {
            throw new IllegalStateException("Unexpected actionable present in fusion evidence: " + actionable);
        }
    }
}
