package com.hartwig.hmftools.serve.sources.vicc;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
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

final class ActionableEvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEvidenceFactory.class);

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet();
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet();

    static {
        RESPONSIVE_DIRECTIONS.add("Responsive");
        RESPONSIVE_DIRECTIONS.add("Sensitivity");
        RESPONSIVE_DIRECTIONS.add("Sensitive");

        RESISTANT_DIRECTIONS.add("Resistant");

        DIRECTIONS_TO_IGNORE.add("Positive");
        DIRECTIONS_TO_IGNORE.add("Negative");
        DIRECTIONS_TO_IGNORE.add("Adverse response");
        DIRECTIONS_TO_IGNORE.add("No responsive");
        DIRECTIONS_TO_IGNORE.add("Not applicable");
        DIRECTIONS_TO_IGNORE.add("Conflicting");
        DIRECTIONS_TO_IGNORE.add("Na");
        DIRECTIONS_TO_IGNORE.add("N/a");
        DIRECTIONS_TO_IGNORE.add("Uncertain significance");
        DIRECTIONS_TO_IGNORE.add("Pathogenic");
        DIRECTIONS_TO_IGNORE.add("Likely pathogenic");
        DIRECTIONS_TO_IGNORE.add("Better outcome");
        DIRECTIONS_TO_IGNORE.add("No benefit");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (myelosupression)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (ototoxicity)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (hyperbilirubinemia)");
        DIRECTIONS_TO_IGNORE.add("Increased toxicity (haemolytic anemia)");
        DIRECTIONS_TO_IGNORE.add("Unknown");
    }

    private ActionableEvidenceFactory() {
    }

    @Nullable
    public static ActionableEvent toActionableEvent(@NotNull ViccEntry entry) {
        String treatment = reformatDrugLabels(entry.association().drugLabels());

        String cancerType = null;
        String doid = null;
        Phenotype phenotype = entry.association().phenotype();
        if (phenotype != null) {
            cancerType = phenotype.description();
            PhenotypeType type = phenotype.type();
            if (type != null) {
                doid = extractDoid(type.id());
            }
        }
        EvidenceLevel level = resolveLevel(entry.association().evidenceLabel());
        EvidenceDirection direction = resolveDirection(entry.association().responseType());

        String url = null;
        EvidenceInfo info = entry.association().evidence().info();
        if (info != null && !info.publications().isEmpty()) {
            url = info.publications().get(0);
        }

        // Consider CancerType, DOID and URL to be optional.
        //TODO: what todo when direction is not correct filled in?
        if (treatment != null && level != null) {
            return ImmutableActionableEvidence.builder()
                    .source(fromViccSource(entry.source()))
                    .treatment(treatment)
                    .cancerType(nullToEmpty(cancerType))
                    .doid(nullToEmpty(doid))
                    .level(level)
                    .direction(nullToEmpty(direction))
                    .url(nullToEmpty(url))
                    .build();
        } else {
            return null;
        }
    }

    @NotNull
    private static EvidenceDirection nullToEmpty(@Nullable EvidenceDirection evidenceDirection) {
        return evidenceDirection != null ? evidenceDirection : EvidenceDirection.NA;
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

    @NotNull
    private static String nullToEmpty(@Nullable String string) {
        return string != null ? string : Strings.EMPTY;
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
                return Strings.EMPTY;
            }
        } else {
            LOGGER.warn("Unexpected DOID string: '{}'", doidString);
            return Strings.EMPTY;
        }
    }

    @NotNull
    private static Knowledgebase fromViccSource(@NotNull ViccSource source) {
        switch (source) {
            case CIVIC:
                return Knowledgebase.CIVIC;
            case CGI:
                return Knowledgebase.CGI;
            case JAX:
                return Knowledgebase.JAX;
            case ONCOKB:
                return Knowledgebase.ONCOKB;
            default:
                throw new IllegalStateException("Source not supported by SERVE: " + source);
        }
    }
}
