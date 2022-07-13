package com.hartwig.hmftools.serve.treatementapproach.filter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class FilterRelevantTreatmentApproachFilter {

    private static final String FIELD_DELIMITER = "\t";

    @NotNull
    public static List<FilterRelevantTreatmentApproachEntry> read(@NotNull String drugClassCurationTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(drugClassCurationTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static List<FilterRelevantTreatmentApproachEntry> fromLines(@NotNull List<String> lines) {
        List<FilterRelevantTreatmentApproachEntry> curatedEntries = Lists.newArrayList();

        for (String line : lines) {
            String[] values = line.split(FIELD_DELIMITER);
            curatedEntries.add(ImmutableFilterRelevantTreatmentApproachEntry.builder()
                    .treatment(values[0])
                    .eventMatch(values[1])
                    .direction(EvidenceDirection.valueOf(values[2]))
                    .level(EvidenceLevel.fromString(values[3]))
                    .build());
        }
        return curatedEntries;
    }
}