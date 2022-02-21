package com.hartwig.hmftools.serve.sources.ckb;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.blacklisting.ImmutableTumorLocationBlacklisting;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklist;
import com.hartwig.hmftools.serve.blacklisting.TumorLocationBlacklisting;

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
    @VisibleForTesting
    public static Set<ActionableEntry> toActionableEntries(@NotNull CkbEntry entry, @NotNull String rawInput) {
        Set<ActionableEntry> actionableEntries = Sets.newHashSet();

        for (Evidence evidence : entry.evidences()) {
            if (hasUsableEvidenceType(evidence.evidenceType())) {
                EvidenceLevel level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
                EvidenceDirection direction = resolveDirection(evidence.responseType());
                String doid = extractDoid(evidence.indication().termId());

                if (level != null && direction != null && doid != null) {
                    String treatment = evidence.therapy().therapyName();
                    String cancerType = evidence.indication().name();

                    Set<String> urls = Sets.newHashSet();
                    for (Reference reference : evidence.references()) {
                        if (reference.url() != null) {
                            urls.add(reference.url());
                        }
                    }

                    int molecularProfileId = entry.profileId();
                    String doidKb = extractDoidKB(evidence.indication().termId());

                    String responseType = extractResponseType(evidence.responseType());

                    Set<String> sourceLinks = Sets.newHashSet();
                    for (Drug drug : evidence.therapy().drugs()) {
                        sourceLinks.add(
                                "https://ckbhome.jax.org/profileResponse/advancedEvidenceFind?molecularProfileId=" + molecularProfileId
                                        + "&drugId=" + drug.id() + "&doId=" + doidKb + "&responseType=" + responseType + "&evidenceType="
                                        + evidence.evidenceType());
                    }

                    Set<TumorLocationBlacklisting> tumorLocationBlacklistings = Sets.newHashSet();
                    tumorLocationBlacklistings.add(ImmutableTumorLocationBlacklisting.builder()
                            .blacklistCancerType(doid.equals("162") ? "Hematologic cancer" : Strings.EMPTY)
                            .blacklistedDoid(doid.equals("162") ? "2531" : Strings.EMPTY)
                            .build());
                    String tumorLocationBlacklist = TumorLocationBlacklist.extractTumorLocationBlacklisting(tumorLocationBlacklistings);
                    String tumorLocationBlacklistDoid = TumorLocationBlacklist.extractTumorLocationDoid(tumorLocationBlacklistings);

                    actionableEntries.add(ImmutableActionableEntry.builder()
                            .sourceEvent(rawInput)
                            .source(Knowledgebase.CKB)
                            .treatment(treatment)
                            .cancerType(cancerType)
                            .doid(doid)
                            .blacklistCancerType(tumorLocationBlacklist)
                            .blacklistedDoid(tumorLocationBlacklistDoid)
                            .level(level)
                            .direction(direction)
                            .sourceUrls(sourceLinks)
                            .evidenceUrls(urls)
                            .build());
                }
            }
        }
        return actionableEntries;
    }

    @NotNull
    @VisibleForTesting
    public static String extractResponseType(@NotNull String responseType) {
        if (responseType.equals("predicted - sensitive")) {
            return "predicted+-+sensitive";
        } else if (responseType.equals("predicted - resistant")) {
            return "predicted+-+resistant";
        } else {
            return responseType;
        }
    }

    @Nullable
    @VisibleForTesting
    static String extractDoidKB(@Nullable String doidString) {
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
    static String extractDoid(@Nullable String doidString) {
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
                if (id.equals("10000003")) {
                    // CKB uses this as Advanced Solid Tumor
                    return "162";
                } else if (id.equals("10000009")) {
                    // CKB uses this as Squamous Cell Carcinoma of Unknown Primary
                    return "1749";
                } else if (id.equals("10000008")) {
                    // CKB uses this as Adenocarcinoma of Unknown Primary
                    return "299";
                } else {
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
    public static boolean hasUsableEvidenceType(@NotNull String evidenceType) {
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
    public static EvidenceLevel resolveLevel(@Nullable String evidenceLabel) {
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
