package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.sources.ckb.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.ckb.curation.EvidenceLevelCurator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ActionableEvidenceFactory {
    private static final Logger LOGGER = LogManager.getLogger(ActionableEvidenceFactory.class);

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet();
    private static final Set<String> EVIDENCE_TYPE = Sets.newHashSet();
    private static final Set<String> EVIDENCE_TYPE_TO_IGNORE = Sets.newHashSet();

    static {
        RESPONSIVE_DIRECTIONS.add("sensitive");
        RESPONSIVE_DIRECTIONS.add("predicted - sensitive");
        RESPONSIVE_DIRECTIONS.add("decreased response");

        RESISTANT_DIRECTIONS.add("resistant");
        RESISTANT_DIRECTIONS.add("predicted - resistant");

        DIRECTIONS_TO_IGNORE.add("unknown");
        DIRECTIONS_TO_IGNORE.add("not applicable");
        DIRECTIONS_TO_IGNORE.add("conflicting");
        DIRECTIONS_TO_IGNORE.add("no benefit");
        DIRECTIONS_TO_IGNORE.add("not predictive");

        EVIDENCE_TYPE.add("Actionable");

        EVIDENCE_TYPE_TO_IGNORE.add("Prognostic");
        EVIDENCE_TYPE_TO_IGNORE.add("Emerging");
        EVIDENCE_TYPE_TO_IGNORE.add("Risk Factor");
        EVIDENCE_TYPE_TO_IGNORE.add("Diagnostic");

    }

    @NotNull
    private final DoidLookup missingDoidLookup;
    @NotNull
    private final DrugCurator drugCurator;
    @NotNull
    private final EvidenceLevelCurator evidenceLevelCurator;

    ActionableEvidenceFactory(@NotNull final DoidLookup missingDoidLookup, @NotNull final DrugCurator drugCurator,
            @NotNull final EvidenceLevelCurator evidenceLevelCurator) {
        this.missingDoidLookup = missingDoidLookup;
        this.drugCurator = drugCurator;
        this.evidenceLevelCurator = evidenceLevelCurator;
    }

    @NotNull
    public Set<ActionableEvent> toActionableEvents(@NotNull CkbEntry entry) {
        Set<ActionableEvent> actionableEvents = Sets.newHashSet();

        for (Evidence evidence : entry.evidences()) {

            if (resolveEvidenceType(evidence.evidenceType())) {
                String therapyName = evidence.therapy().therapyName();

                EvidenceLevel level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
                EvidenceDirection direction = resolveDirection(evidence.responseType());
                String cancerType = evidence.indication().name();
                String doid = evidence.indication().termId();
                Set<String> urls = Sets.newHashSet();

                Map<String, Set<String>> cancerTypeToDoidsMap = buildCancerTypeToDoidsMap(cancerType, extractDoid(doid));

                for (Reference reference : evidence.references()) {
                    if (reference.url() != null) {
                        urls.add(reference.url());
                    }
                }

                // TODO: Determine if curation of drugs and evidence is needed

                if (therapyName != null && level != null && direction != null) {
                    ImmutableActionableEvidence.Builder builder =
                            ImmutableActionableEvidence.builder().source(Knowledgebase.CKB).level(level).direction(direction).urls(urls);
                    for (Map.Entry<String, Set<String>> cancerTypeEntry : cancerTypeToDoidsMap.entrySet()) {
                        for (String doidSet : cancerTypeEntry.getValue()) {
                            actionableEvents.add(builder.cancerType(cancerType).doid(doidSet).treatment(therapyName).build());
                        }
                    }
                }
            }
        }

        return actionableEvents;
    }

    @Nullable
    @VisibleForTesting
    static String extractDoid(@Nullable String doidString) {
        if (doidString == null) {
            return null;
        }

        String[] parts = doidString.split(":");
        if (parts.length == 2) {
            if (parts[0].equalsIgnoreCase("doid")) {
                return parts[1];
            } else if (parts[0].equalsIgnoreCase("jax")) {
                return doidString;
            } else {
                return null;
            }
        } else {
            LOGGER.warn("Unexpected Doid string: '{}'", doidString);
            return null;
        }
    }

    @NotNull
    private Set<String> lookupDoids(@NotNull String cancerType) {
        Set<String> doids = missingDoidLookup.lookupDoidsForCancerType(cancerType);
        if (doids != null) {
            return doids;
        } else {
            LOGGER.warn("Could not resolve doids for VICC cancer type '{}'", cancerType);
            return Sets.newHashSet();
        }
    }

    @NotNull
    private Map<String, Set<String>> buildCancerTypeToDoidsMap(@Nullable String cancerType, @Nullable String doid) {
        Map<String, Set<String>> cancerTypeToDoidsMap = Maps.newHashMap();
        if (cancerType != null) {
            if (doid != null) {
                if (doid.split(";")[0].equalsIgnoreCase("jax")) {
                    cancerTypeToDoidsMap.put(cancerType, Sets.newHashSet("0050686"));
                } else {
                    cancerTypeToDoidsMap.put(cancerType, Sets.newHashSet(doid));

                }
            }
        }
        return cancerTypeToDoidsMap;
    }

    public void evaluateCuration() {
        drugCurator.reportUnusedCurationKeys();
        evidenceLevelCurator.reportUnusedCurationKeys();
    }

    private static boolean resolveEvidenceType(@NotNull String evidenceType) {
        boolean usableEvidence = false;
        if (EVIDENCE_TYPE.contains(evidenceType)) {
            usableEvidence = true;
        } else if (EVIDENCE_TYPE_TO_IGNORE.contains(evidenceType)) {
            usableEvidence = false;
        }
        return usableEvidence;
    }

    @Nullable
    @$VisibleForTesting
    static EvidenceLevel resolveLevel(@Nullable String evidenceLabel) {
        if (evidenceLabel == null) {
            return null;
        }

        EvidenceLevel level = EvidenceLevel.fromString(evidenceLabel);
        if (level == null) {
            LOGGER.warn("Could not resolve evidence label '{}'", evidenceLabel);
        }
        return level;
    }

    @Nullable
    @VisibleForTesting
    static EvidenceDirection resolveDirection(@Nullable String direction) {
        if (direction == null) {
            return null;
        }

        if (RESPONSIVE_DIRECTIONS.contains(direction)) {
            return EvidenceDirection.RESPONSIVE;
        } else if (RESISTANT_DIRECTIONS.contains(direction)) {
            return EvidenceDirection.RESISTANT;
        }

        if (!DIRECTIONS_TO_IGNORE.contains(direction)) {
            LOGGER.warn("Could not resolve CKB direction '{}'", direction);
        }
        return null;
    }
}
