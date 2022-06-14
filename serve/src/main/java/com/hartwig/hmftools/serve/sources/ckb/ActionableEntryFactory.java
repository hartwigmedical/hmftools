package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproaches.RelevantTreatmentApproaches;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.cancertype.CancerType;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.treatment.ImmutableTreatment;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ActionableEntryFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEntryFactory.class);

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> PREDICTED_RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> PREDICTED_RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet();

    private static final Set<String> USABLE_EVIDENCE_TYPES = Sets.newHashSet();
    private static final Set<String> EVIDENCE_TYPES_TO_IGNORE = Sets.newHashSet();

    static {
        RESPONSIVE_DIRECTIONS.add("sensitive");
        PREDICTED_RESPONSIVE_DIRECTIONS.add("predicted - sensitive");

        RESISTANT_DIRECTIONS.add("resistant");
        PREDICTED_RESISTANT_DIRECTIONS.add("predicted - resistant");

        DIRECTIONS_TO_IGNORE.add("unknown");
        DIRECTIONS_TO_IGNORE.add("not applicable");
        DIRECTIONS_TO_IGNORE.add("conflicting");
        DIRECTIONS_TO_IGNORE.add("no benefit");
        DIRECTIONS_TO_IGNORE.add("not predictive");
        DIRECTIONS_TO_IGNORE.add("decreased response");

        USABLE_EVIDENCE_TYPES.add("Actionable");

        EVIDENCE_TYPES_TO_IGNORE.add("Prognostic");
        EVIDENCE_TYPES_TO_IGNORE.add("Emerging");
        EVIDENCE_TYPES_TO_IGNORE.add("Risk Factor");
        EVIDENCE_TYPES_TO_IGNORE.add("Diagnostic");
    }

    ActionableEntryFactory() {
    }

    @NotNull
    public static Set<ActionableEntry> toActionableEntries(@NotNull CkbEntry entry, @NotNull String sourceEvent) {
        Set<ActionableEntry> actionableEntries = Sets.newHashSet();
        Set<String> drugClasses = Sets.newHashSet();
        for (Evidence evidence : evidencesWithUsableType(entry.evidences())) {
            EvidenceLevel level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
            EvidenceDirection direction = resolveDirection(evidence.responseType());
            String doid = extractAndCurateDoid(evidence.indication().termId());

            for (RelevantTreatmentApproaches relevantTreatmentApproaches : evidence.relevantTreatmentApproaches()) {
                DrugClass drugClass = relevantTreatmentApproaches.drugClass();
                drugClasses.add(drugClass == null ? Strings.EMPTY: drugClass.drugClass());
            }

            if (level != null && direction != null && doid != null) {
                String treatment = evidence.therapy().therapyName();
                String cancerType = evidence.indication().name();

                Set<String> evidenceUrls = Sets.newHashSet();
                for (Reference reference : evidence.references()) {
                    if (reference.url() != null) {
                        evidenceUrls.add(reference.url());
                    }
                }

                String sourceCancerTypeId = extractSourceCancerTypeId(evidence.indication().termId());

                String responseType = toUrlString(evidence.responseType());

                Set<String> sourceUrls = Sets.newHashSet();
                for (Drug drug : evidence.therapy().drugs()) {
                    sourceUrls.add("https://ckbhome.jax.org/profileResponse/advancedEvidenceFind?molecularProfileId=" + entry.profileId()
                            + "&drugId=" + drug.id() + "&doId=" + sourceCancerTypeId + "&responseType=" + responseType + "&evidenceType="
                            + evidence.evidenceType());
                }

                Set<CancerType> blacklistedCancerTypes = Sets.newHashSet();
                if (doid.equals("162")) {
                    blacklistedCancerTypes.add(ImmutableCancerType.builder().name("Leukemia").doid("1240").build());
                    blacklistedCancerTypes.add(ImmutableCancerType.builder().name("Refractory hematologic cancer").doid("712").build());
                    blacklistedCancerTypes.add(ImmutableCancerType.builder().name("Bone marrow cancer").doid("4960").build());
                }

                actionableEntries.add(ImmutableActionableEntry.builder()
                        .source(Knowledgebase.CKB)
                        .sourceEvent(sourceEvent)
                        .sourceUrls(sourceUrls)
                        .treatment(ImmutableTreatment.builder().treament(treatment).drugClasses(drugClasses).build())
                        .applicableCancerType(ImmutableCancerType.builder().name(cancerType).doid(doid).build())
                        .blacklistCancerTypes(blacklistedCancerTypes)
                        .level(level)
                        .direction(direction)
                        .evidenceUrls(evidenceUrls)
                        .build());
            }
        }
        return actionableEntries;
    }

    @NotNull
    private static List<Evidence> evidencesWithUsableType(@NotNull List<Evidence> evidences) {
        List<Evidence> filtered = Lists.newArrayList();
        for (Evidence evidence : evidences) {
            if (hasUsableEvidenceType(evidence.evidenceType())) {
                filtered.add(evidence);
            }
        }
        return filtered;
    }

    @NotNull
    @VisibleForTesting
    static String toUrlString(@NotNull String string) {
        return string.replaceAll(" ", "+");
    }

    @Nullable
    @VisibleForTesting
    static String extractSourceCancerTypeId(@Nullable String doidString) {
        if (doidString == null) {
            return null;
        }

        String[] parts = doidString.split(":");
        if (parts.length == 2) {
            String source = parts[0];
            String id = parts[1];
            if (source.equalsIgnoreCase("doid") || source.equalsIgnoreCase("jax")) {
                return id;
            } else {
                return null;
            }
        } else {
            LOGGER.warn("Unexpected DOID string in CKB: '{}'", doidString);
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    static String extractAndCurateDoid(@Nullable String doidString) {
        if (doidString == null) {
            return null;
        }

        String[] parts = doidString.split(":");
        if (parts.length == 2) {
            String source = parts[0];
            String id = parts[1];
            if (source.equalsIgnoreCase("doid")) {
                return id;
            } else if (source.equalsIgnoreCase("jax")) {
                switch (id) {
                    case "10000003":
                        // CKB uses this as Advanced Solid Tumor
                        return "162";
                    case "10000009":
                        // CKB uses this as Squamous Cell Carcinoma of Unknown Primary
                        return "1749";
                    case "10000008":
                        // CKB uses this as Adenocarcinoma of Unknown Primary
                        return "299";
                    default:
                        // CKB uses 10000005 for configuring "Not a cancer". We can ignore these.
                        if (!id.equals("10000005")) {
                            LOGGER.warn("Unexpected DOID string annotated by CKB: '{}'", doidString);
                        }
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

    @VisibleForTesting
    static boolean hasUsableEvidenceType(@NotNull String evidenceType) {
        if (USABLE_EVIDENCE_TYPES.contains(evidenceType)) {
            return true;
        } else {
            if (!EVIDENCE_TYPES_TO_IGNORE.contains(evidenceType)) {
                LOGGER.warn("Unrecognized CKB evidence type: '{}'", evidenceType);
            }
            return false;
        }
    }

    @Nullable
    @VisibleForTesting
    static EvidenceLevel resolveLevel(@Nullable String evidenceLabel) {
        if (evidenceLabel == null || evidenceLabel.equals("NA")) {
            return null;
        }

        EvidenceLevel level = EvidenceLevel.fromString(evidenceLabel);
        if (level == null) {
            LOGGER.warn("Could not resolve CKB evidence level: '{}'", evidenceLabel);
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
        } else if (PREDICTED_RESPONSIVE_DIRECTIONS.contains(direction)) {
            return EvidenceDirection.PREDICTED_RESPONSIVE;
        } else if (RESISTANT_DIRECTIONS.contains(direction)) {
            return EvidenceDirection.RESISTANT;
        } else if (PREDICTED_RESISTANT_DIRECTIONS.contains(direction)) {
            return EvidenceDirection.PREDICTED_RESISTANT;
        }

        if (!DIRECTIONS_TO_IGNORE.contains(direction)) {
            LOGGER.warn("Could not resolve CKB direction '{}'", direction);
        }
        return null;
    }
}