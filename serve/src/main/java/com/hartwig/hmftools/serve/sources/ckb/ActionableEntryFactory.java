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
import com.hartwig.hmftools.common.serve.actionability.ImmutableTreatment;
import com.hartwig.hmftools.common.serve.cancertype.CancerType;
import com.hartwig.hmftools.common.serve.cancertype.CancerTypeConstants;
import com.hartwig.hmftools.common.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.ckb.treatementapproach.ImmutableRelevantTreatmentApprochCurationEntryKey;
import com.hartwig.hmftools.serve.sources.ckb.treatementapproach.RelevantTreatmentApprochCurationEntryKey;
import com.hartwig.hmftools.serve.sources.ckb.treatementapproach.RelevantTreatmentAproachCuration;

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
    public static Set<ActionableEntry> toActionableEntries(@NotNull CkbEntry entry, @NotNull String sourceEvent,
            @NotNull RelevantTreatmentAproachCuration curator, @NotNull String gene, @NotNull EventType eventType) {
        Set<ActionableEntry> actionableEntries = Sets.newHashSet();

        for (Evidence evidence : evidencesWithUsableType(entry.evidences())) {
            EvidenceLevel level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
            EvidenceDirection direction = resolveDirection(evidence.responseType());
            String[] sourceCancerTypes = extractSourceCancerTypeId(evidence.indication().termId());
            String doid = extractAndCurateDoid(sourceCancerTypes);

            if (level != null && direction != null && doid != null) {
                String treatment = evidence.therapy().therapyName();
                String cancerType = evidence.indication().name();

                Set<String> evidenceUrls = Sets.newHashSet();
                for (Reference reference : evidence.references()) {
                    if (reference.url() != null) {
                        evidenceUrls.add(reference.url());
                    }
                }

                String sourceCancerTypeId;
                if (extractSourceCancerTypeId(evidence.indication().termId()) == null) {
                    sourceCancerTypeId = null;
                } else {
                    assert extractSourceCancerTypeId(evidence.indication().termId()) != null;
                    assert extractSourceCancerTypeId(evidence.indication().termId()).length == 2;
                    sourceCancerTypeId = extractSourceCancerTypeId(evidence.indication().termId())[1];
                }

                String responseType = toUrlString(evidence.responseType());

                Set<String> sourceUrls = Sets.newHashSet();
                for (Drug drug : evidence.therapy().drugs()) {
                    sourceUrls.add("https://ckbhome.jax.org/profileResponse/advancedEvidenceFind?molecularProfileId=" + entry.profileId()
                            + "&drugId=" + drug.id() + "&doId=" + sourceCancerTypeId + "&responseType=" + responseType + "&evidenceType="
                            + evidence.evidenceType());
                }

                Set<CancerType> blacklistedCancerTypes = Sets.newHashSet();
                if (doid.equals(CancerTypeConstants.CANCER_DOID)) {
                    blacklistedCancerTypes.add(CancerTypeConstants.LEUKEMIA_TYPE);
                    blacklistedCancerTypes.add(CancerTypeConstants.REFRACTORY_HEMATOLOGIC_TYPE);
                    blacklistedCancerTypes.add(CancerTypeConstants.BONE_MARROW_TYPE);
                }

                Set<String> sourceRelevantTreatmentApproaches = Sets.newHashSet();
                Set<String> curatedRelevantTreatmentApproaches = Sets.newHashSet();
                for (RelevantTreatmentApproaches relevantTreatmentApproaches : evidence.relevantTreatmentApproaches()) {
                    DrugClass relevantTreatmentApproachesInfo = relevantTreatmentApproaches.drugClass();

                    if (relevantTreatmentApproachesInfo != null) {
                        sourceRelevantTreatmentApproaches.add(relevantTreatmentApproachesInfo.drugClass());
                    }
                }

                String treatmentApproachString = String.join(",", sourceRelevantTreatmentApproaches);
                String treatmentApproachInterpret = Strings.EMPTY;
                if (sourceRelevantTreatmentApproaches.isEmpty()) {
                    treatmentApproachInterpret = null;
                } else if (treatmentApproachString.substring(treatmentApproachString.length() - 1).equals(",")) {
                    treatmentApproachInterpret = treatmentApproachString.substring(0, treatmentApproachString.length() - 1);
                } else {
                    treatmentApproachInterpret = treatmentApproachString;
                }

                RelevantTreatmentApprochCurationEntryKey key = ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                        .treatment(treatment)
                        .treatmentApproach(treatmentApproachInterpret == null || treatmentApproachInterpret.isEmpty()
                                ? null
                                : treatmentApproachInterpret)
                        .event(gene + " " + eventType)
                        .direction(direction)
                        .build();

                curatedRelevantTreatmentApproaches.add(curator.isMatch(key));

                actionableEntries.add(ImmutableActionableEntry.builder()
                        .source(Knowledgebase.CKB)
                        .sourceEvent(sourceEvent)
                        .sourceUrls(sourceUrls)
                        .treatment(ImmutableTreatment.builder()
                                .treament(treatment)
                                .sourceRelevantTreatmentApproaches(sourceRelevantTreatmentApproaches)
                                .relevantTreatmentApproaches(curatedRelevantTreatmentApproaches)
                                .build())
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
    static String[] extractSourceCancerTypeId(@Nullable String doidString) {
        if (doidString == null) {
            return null;
        }

        String[] parts = doidString.split(":");
        if (parts.length == 2) {
            String source = parts[0];
            if (source.equalsIgnoreCase("doid") || source.equalsIgnoreCase("jax")) {
                return parts;
            } else {
                LOGGER.warn("Unexpected length of doid string '{}'", doidString);
                return null;
            }
        } else {
            LOGGER.warn("Unexpected DOID string in CKB: '{}'", doidString);
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    static String extractAndCurateDoid(@Nullable String[] doidString) {
        if (doidString == null) {
            return null;
        }

        assert doidString.length == 2;
        String source = doidString[0];
        String id = doidString[1];
        if (source.equalsIgnoreCase("doid")) {
            return id;
        } else if (source.equalsIgnoreCase("jax")) {
            switch (id) {
                case CancerTypeConstants.JAX_ADVANCES_SOLID_TUMORS:
                    // CKB uses this as Advanced Solid Tumor
                    return CancerTypeConstants.CANCER_DOID;
                case CancerTypeConstants.JAX_SQUAMOUD_CELL_CARCINOMA_OF_UNKNOWN_PRIMARY:
                    // CKB uses this as Squamous Cell Carcinoma of Unknown Primary
                    return CancerTypeConstants.SQUAMOUD_CELL_CARCINOMA_OF_UNKNOWN_PRIMARY;
                case CancerTypeConstants.JAX_ADENOCARCINOMA_OF_UNKNOWN_PRIMARY:
                    // CKB uses this as Adenocarcinoma of Unknown Primary
                    return CancerTypeConstants.ADENOCARCINOMA_OF_UNKNOWN_PRIMARY;
                default:
                    // CKB uses 10000005 for configuring "Not a cancer". We can ignore these.
                    if (!id.equals(CancerTypeConstants.JAX_NOT_CANCER)) {
                        LOGGER.warn("Unexpected DOID string annotated by CKB: '{}'", source + ":" + id);
                    }
                    return null;
            }
        } else {
            LOGGER.warn("Unexpected source '{}'", source);
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