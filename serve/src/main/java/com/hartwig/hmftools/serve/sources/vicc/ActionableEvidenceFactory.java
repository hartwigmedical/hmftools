package com.hartwig.hmftools.serve.sources.vicc;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ActionableEvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEvidenceFactory.class);

    private static final String CANCER_TYPE_SEPARATOR = ";";

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet();

    static {
        RESPONSIVE_DIRECTIONS.add("Responsive");
        RESPONSIVE_DIRECTIONS.add("Sensitivity");
        RESPONSIVE_DIRECTIONS.add("Sensitive");

        RESISTANT_DIRECTIONS.add("Resistant");

        DIRECTIONS_TO_IGNORE.add("Adverse response");
        DIRECTIONS_TO_IGNORE.add("No responsive");
        DIRECTIONS_TO_IGNORE.add("Not applicable");
        DIRECTIONS_TO_IGNORE.add("Conflicting");
        DIRECTIONS_TO_IGNORE.add("Na");
        DIRECTIONS_TO_IGNORE.add("N/a");
        DIRECTIONS_TO_IGNORE.add("No benefit");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (myelosupression)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (ototoxicity)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (hyperbilirubinemia)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (haemolytic anemia)");
        DIRECTIONS_TO_IGNORE.add("Unknown");
    }

    @NotNull
    private final DoidLookup missingDoidLookup;

    public ActionableEvidenceFactory(@NotNull final DoidLookup missingDoidLookup) {
        this.missingDoidLookup = missingDoidLookup;
    }

    @NotNull
    public Set<ActionableEvent> toActionableEvents(@NotNull ViccEntry entry) {
        Set<ActionableEvent> actionableEvents = Sets.newHashSet();

        String treatment = reformatDrugLabels(entry.association().drugLabels());
        EvidenceLevel level = resolveLevel(entry.association().evidenceLabel());

        if (treatment != null && level != null) {
            EvidenceDirection direction = resolveDirection(entry.association().responseType());
            Set<String> urls = resolveUrls(entry.association().evidence().info());

            ImmutableActionableEvidence.Builder builder = ImmutableActionableEvidence.builder()
                    .source(fromViccSource(entry.source()))
                    .treatment(treatment)
                    .level(level)
                    .direction(nullToEmpty(direction))
                    .urls(urls);

            String cancerType = resolveCancerType(entry.association().phenotype());

            if (cancerType != null) {
                if (cancerType.contains(CANCER_TYPE_SEPARATOR)) {
                    String[] parts = cancerType.split(CANCER_TYPE_SEPARATOR);
                    for (String part : parts) {
                        // We always look up the DOIDs when there is aggregate cancer type information as the DOID in this case is unreliable.
                        for (String doid : lookupDoids(part)) {
                            actionableEvents.add(builder.cancerType(part).doid(doid).build());
                        }
                    }
                } else {
                    String doidEntry = resolveDoid(entry.association().phenotype());
                    Set<String> doids;
                    if (doidEntry == null) {
                        doids = lookupDoids(cancerType);
                    } else {
                        doids = Sets.newHashSet(doidEntry);
                    }

                    for (String doid : doids) {
                        actionableEvents.add(builder.cancerType(cancerType).doid(doid).build());
                    }
                }
            }
        }

        return actionableEvents;
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

    @Nullable
    @VisibleForTesting
    static String reformatDrugLabels(@Nullable String drugLabels) {
        if (drugLabels == null) {
            return null;
        }

        String drugSeparator = ",";
        String[] parts = drugLabels.split(drugSeparator);
        StringJoiner joiner = new StringJoiner(drugSeparator);
        for (String part : parts) {
            joiner.add(reformatField(part));
        }
        return joiner.toString();
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
    private static String resolveCancerType(@Nullable Phenotype phenotype) {
        return phenotype != null ? phenotype.description() : null;
    }

    @Nullable
    private static String resolveDoid(@Nullable Phenotype phenotype) {
        if (phenotype != null) {
            PhenotypeType type = phenotype.type();
            if (type != null) {
                return extractDoid(type.id());
            }
        }

        return null;
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
            } else {
                return null;
            }
        } else {
            LOGGER.warn("Unexpected Doid string: '{}'", doidString);
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    static EvidenceDirection resolveDirection(@Nullable String direction) {
        String effectiveDirection = reformatField(direction);
        if (effectiveDirection == null) {
            return null;
        }

        if (RESPONSIVE_DIRECTIONS.contains(effectiveDirection)) {
            return EvidenceDirection.RESPONSIVE;
        } else if (RESISTANT_DIRECTIONS.contains(effectiveDirection)) {
            return EvidenceDirection.RESISTANT;
        }

        if (!DIRECTIONS_TO_IGNORE.contains(effectiveDirection)) {
            LOGGER.warn("Could not resolve VICC direction '{}'", effectiveDirection);
        }
        return null;
    }

    @NotNull
    private static Set<String> resolveUrls(@Nullable EvidenceInfo info) {
        return info != null ? Sets.newHashSet(info.publications()) : Sets.newHashSet();
    }

    @NotNull
    private static EvidenceDirection nullToEmpty(@Nullable EvidenceDirection evidenceDirection) {
        return evidenceDirection != null ? evidenceDirection : EvidenceDirection.NA;
    }

    @NotNull
    private static String nullToEmpty(@Nullable String string) {
        return string != null ? string : Strings.EMPTY;
    }

    @Nullable
    @VisibleForTesting
    static String reformatField(@Nullable String direction) {
        if (direction == null) {
            return null;
        } else if (direction.length() < 2) {
            return direction.toUpperCase();
        } else {
            return direction.substring(0, 1).toUpperCase() + direction.substring(1).toLowerCase();
        }
    }

    @NotNull
    private static Knowledgebase fromViccSource(@NotNull ViccSource source) {
        switch (source) {
            case CIVIC:
                return Knowledgebase.VICC_CIVIC;
            case CGI:
                return Knowledgebase.VICC_CGI;
            case JAX:
                return Knowledgebase.VICC_JAX;
            case ONCOKB:
                return Knowledgebase.VICC_ONCOKB;
            default:
                throw new IllegalStateException("Source not supported by SERVE: " + source);
        }
    }
}
