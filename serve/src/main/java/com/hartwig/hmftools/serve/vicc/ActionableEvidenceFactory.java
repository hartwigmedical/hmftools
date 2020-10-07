package com.hartwig.hmftools.serve.vicc;

import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class ActionableEvidenceFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableEvidenceFactory.class);

    private static final Set<String> RESPONSIVE_DIRECTIONS = Sets.newHashSet("Responsive", "Sensitivity", "Sensitive");
    private static final Set<String> RESISTANT_DIRECTIONS = Sets.newHashSet("Resistant");
    private static final Set<String> DIRECTIONS_TO_IGNORE = Sets.newHashSet("Not applicable", "Conflicting", "No benefit", "NA");

    private ActionableEvidenceFactory() {
    }

    @Nullable
    public static ActionableEvidence toActionableEvidence(@NotNull ViccEntry entry) {
        String drugs = reformatDrugs(entry.association().drugLabels());

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
        String level = entry.association().evidenceLabel();
        EvidenceDirection direction = resolveDirection(entry.association().responseType());

        String url = null;
        EvidenceInfo info = entry.association().evidence().info();
        if (info != null && !info.publications().isEmpty()) {
            url = info.publications().get(0);
        }

        // Consider CancerType, DOID and URL to be optional.
        if (drugs != null && level != null && direction != null) {
            return ImmutableActionableEvidence.builder()
                    .drugs(drugs)
                    .cancerType(nullToEmpty(cancerType))
                    .doid(nullToEmpty(doid))
                    .level(level)
                    .direction(direction)
                    .url(nullToEmpty(url))
                    .build();
        } else {
            return null;
        }
    }

    @NotNull
    private static String nullToEmpty(@Nullable String string) {
        return string != null ? string : Strings.EMPTY;
    }

    @Nullable
    @VisibleForTesting
    static String reformatDrugs(@Nullable String drugs) {
        if (drugs == null) {
            return null;
        }

        String drugSeparator = ",";
        String[] parts = drugs.split(drugSeparator);
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
}
