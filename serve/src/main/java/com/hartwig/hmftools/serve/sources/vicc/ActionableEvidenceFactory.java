package com.hartwig.hmftools.serve.sources.vicc;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;

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

        // These directions only appear in evidence which lacks either level or drugs
        DIRECTIONS_TO_IGNORE.add("Pathogenic");
        DIRECTIONS_TO_IGNORE.add("Likely pathogenic");
        DIRECTIONS_TO_IGNORE.add("Positive");
        DIRECTIONS_TO_IGNORE.add("Negative");
        DIRECTIONS_TO_IGNORE.add("Uncertain significance");
        DIRECTIONS_TO_IGNORE.add("Better outcome");
    }

    @NotNull
    private final DoidLookup missingDoidLookup;
    @NotNull
    private final DrugCurator drugCurator;
    @NotNull
    private final EvidenceLevelCurator evidenceLevelCurator;

    public ActionableEvidenceFactory(@NotNull final DoidLookup missingDoidLookup, @NotNull final DrugCurator drugCurator,
            @NotNull final EvidenceLevelCurator evidenceLevelCurator) {
        this.missingDoidLookup = missingDoidLookup;
        this.drugCurator = drugCurator;
        this.evidenceLevelCurator = evidenceLevelCurator;
    }

    @NotNull
    public Set<ActionableEvent> toActionableEvents(@NotNull ViccEntry entry, @NotNull String rawInput) {
        Set<ActionableEvent> actionableEvents = Sets.newHashSet();

        boolean isSupportive = isSupportiveEntry(entry);
        String treatment = reformatDrugLabels(entry.association().drugLabels());
        EvidenceLevel level = resolveLevel(entry.association().evidenceLabel());
        EvidenceDirection direction = resolveDirection(entry.association().responseType());

        if (isSupportive && treatment != null && level != null && direction != null) {
            Set<String> urls = resolveUrls(entry.association().evidence().info());

            Map<String, Set<String>> cancerTypeToDoidsMap = buildCancerTypeToDoidsMap(resolveCancerType(entry.association().phenotype()),
                    resolveDoid(entry.association().phenotype()));

            level = evidenceLevelCurator.curate(entry.source(), entry.genes(), treatment, level, direction);
            List<List<String>> drugLists = drugCurator.curate(entry.source(), level, treatment);

            ImmutableActionableEvidence.Builder builder = ImmutableActionableEvidence.builder()
                    .rawInput(rawInput)
                    .source(fromViccSource(entry.source()))
                    .level(level)
                    .direction(direction)
                    .urlSource(Strings.EMPTY)
                    .urls(urls);

            for (Map.Entry<String, Set<String>> cancerTypeEntry : cancerTypeToDoidsMap.entrySet()) {
                String cancerType = cancerTypeEntry.getKey();
                for (String doid : cancerTypeEntry.getValue()) {
                    for (List<String> drugList : drugLists) {
                        actionableEvents.add(builder.cancerType(cancerType).doid(doid).treatment(formatDrugList(drugList)).build());
                    }
                }
            }
        }

        return actionableEvents;
    }

    private static boolean isSupportiveEntry(@NotNull ViccEntry entry) {
        // CIViC contributes entries that seem "sensitive" or "resistant" but are not "supportive" and rather do not support the evidence.
        // We do not want to generate actionability for them (see also INC-92)
        if (entry.kbSpecificObject() instanceof Civic) {
            String direction = ((Civic) entry.kbSpecificObject()).evidenceItem().evidenceDirection();
            if (direction == null || direction.equals("Supports")) {
                return true;
            } else if (direction.equals("Does Not Support")) {
                return false;
            } else {
                LOGGER.warn("Unrecognized CIViC direction entry '{}' in entry {}", direction, entry);
                return true;
            }
        } else {
            return true;
        }
    }

    public void evaluateCuration() {
        drugCurator.reportUnusedCurationKeys();
        evidenceLevelCurator.reportUnusedCurationKeys();
    }

    @NotNull
    private static String formatDrugList(@NotNull List<String> drugList) {
        List<String> sortedDrugs = Lists.newArrayList(drugList);
        sortedDrugs.sort(Comparator.naturalOrder());

        StringJoiner joiner = new StringJoiner(" + ");
        for (String drug : sortedDrugs) {
            joiner.add(drug);
        }
        return joiner.toString();
    }

    @NotNull
    private Map<String, Set<String>> buildCancerTypeToDoidsMap(@Nullable String cancerType, @Nullable String doid) {
        Map<String, Set<String>> cancerTypeToDoidsMap = Maps.newHashMap();
        if (cancerType != null) {
            if (cancerType.contains(CANCER_TYPE_SEPARATOR)) {
                String[] parts = cancerType.split(CANCER_TYPE_SEPARATOR);
                for (String part : parts) {
                    // We always look up the DOIDs when there is aggregate cancer type information as the DOID in this case is unreliable.
                    cancerTypeToDoidsMap.put(part, lookupDoids(part));
                }
            } else if (doid != null) {
                cancerTypeToDoidsMap.put(cancerType, Sets.newHashSet(doid));
            } else {
                cancerTypeToDoidsMap.put(cancerType, lookupDoids(cancerType));
            }
        }
        return cancerTypeToDoidsMap;
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

    @Nullable
    @VisibleForTesting
    static String reformatField(@Nullable String field) {
        if (field == null) {
            return null;
        } else if (field.length() < 2) {
            return field.toUpperCase();
        } else {
            return field.substring(0, 1).toUpperCase() + field.substring(1).toLowerCase();
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
