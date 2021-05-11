package com.hartwig.hmftools.serve.sources.ckb;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ActionableEntryFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEntryFactory.class);

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet();

    private static final Set<String> USABLE_EVIDENCE_TYPES = Sets.newHashSet();
    private static final Set<String> EVIDENCE_TYPES_TO_IGNORE = Sets.newHashSet();

    static {
        RESPONSIVE_DIRECTIONS.add("sensitive");
        RESPONSIVE_DIRECTIONS.add("predicted - sensitive");

        RESISTANT_DIRECTIONS.add("resistant");
        RESISTANT_DIRECTIONS.add("predicted - resistant");

        DIRECTIONS_TO_IGNORE.add("unknown");
        DIRECTIONS_TO_IGNORE.add("not applicable");
        DIRECTIONS_TO_IGNORE.add("conflicting");
        DIRECTIONS_TO_IGNORE.add("no benefit");
        DIRECTIONS_TO_IGNORE.add("not predictive");

        // TODO: Determine first what they mean with below direction
        DIRECTIONS_TO_IGNORE.add("decreased response");

        USABLE_EVIDENCE_TYPES.add("Actionable");

        // TODO: Determine if all those can be ignored
        EVIDENCE_TYPES_TO_IGNORE.add("Prognostic");
        EVIDENCE_TYPES_TO_IGNORE.add("Emerging");
        EVIDENCE_TYPES_TO_IGNORE.add("Risk Factor");
        EVIDENCE_TYPES_TO_IGNORE.add("Diagnostic");
    }

    ActionableEntryFactory() {
    }

    @NotNull
    public Set<ActionableEntry> toActionableEntries(@NotNull CkbEntry entry) {
        Set<ActionableEntry> actionableEntries = Sets.newHashSet();

        for (Evidence evidence : entry.evidences()) {
            if (hasUsableEvidenceType(evidence)) {
                String treatment = evidence.therapy().therapyName();
                EvidenceLevel level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
                EvidenceDirection direction = resolveDirection(evidence.responseType());
                String cancerType = evidence.indication().name();
                String doid = extractDoid(evidence.indication().termId());

                Set<String> urls = Sets.newHashSet();
                for (Reference reference : evidence.references()) {
                    if (reference.url() != null) {
                        urls.add(reference.url());
                    }
                }

                if (level != null && direction != null && doid != null) {
                    actionableEntries.add(ImmutableActionableEntry.builder()
                            .source(Knowledgebase.CKB)
                            .treatment(treatment)
                            .cancerType(cancerType)
                            .doid(doid)
                            .level(level)
                            .direction(direction)
                            .urls(urls)
                            .build());
                }
            }
        }

        return actionableEntries;
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
                if (parts[1].equals("10000003")) {
                    return "0050686";
                } else {
                    LOGGER.warn("Unexpected DOID string annotated by JAX: '{}'", doidString);
                    return null;
                }
            } else {
                return null;
            }
        } else {
            LOGGER.warn("Unexpected DOID string in CKB: '{}'", doidString);
            return null;
        }
    }

    private static boolean hasUsableEvidenceType(@NotNull Evidence evidence) {
        String evidenceType = evidence.evidenceType();
        if (USABLE_EVIDENCE_TYPES.contains(evidenceType)) {
            return true;
        } else {
            if (!EVIDENCE_TYPES_TO_IGNORE.contains(evidenceType)) {
                LOGGER.warn("Unrecognized CKB evidence type '{}'", evidenceType);
            }
            return false;
        }
    }

    @Nullable
    static EvidenceLevel resolveLevel(@Nullable String evidenceLabel) {
        if (evidenceLabel == null || evidenceLabel.equals("NA")) {
            return null;
        }

        EvidenceLevel level = EvidenceLevel.fromString(evidenceLabel);
        if (level == null) {
            LOGGER.warn("Could not resolve CKB evidence level '{}'", evidenceLabel);
        }
        return level;
    }

    @Nullable
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
