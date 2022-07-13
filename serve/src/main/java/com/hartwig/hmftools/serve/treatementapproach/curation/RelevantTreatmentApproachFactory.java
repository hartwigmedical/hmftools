package com.hartwig.hmftools.serve.treatementapproach.curation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class RelevantTreatmentApproachFactory {

    private static final String FIELD_DELIMITER = "\t";

    @NotNull
    public static Map<RelevantTreatmentApproachKey, RelevantTreatmentApproch> read(@NotNull String drugClassCurationTsv)
            throws IOException {
        List<String> lines = Files.readAllLines(new File(drugClassCurationTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static Map<RelevantTreatmentApproachKey, RelevantTreatmentApproch> fromLines(@NotNull List<String> lines) {
        Map<RelevantTreatmentApproachKey, RelevantTreatmentApproch> curatedEntries = Maps.newHashMap();

        for (String line : lines) {
            String[] values = line.split(FIELD_DELIMITER);
            RelevantTreatmentApproachKey drugClassKey = ImmutableRelevantTreatmentApproachKey.builder()
                    .treatment(values[0])
                    .treatmentApproach(values[1])
                    .matchEvent(values[2])
                    .build();
            curatedEntries.put(drugClassKey,
                    ImmutableRelevantTreatmentApproch.builder()
                            .treatmentApproachKey(drugClassKey)
                            .curatedtreatmentApproach(values[3])
                            .build());
        }
        return curatedEntries;
    }
}