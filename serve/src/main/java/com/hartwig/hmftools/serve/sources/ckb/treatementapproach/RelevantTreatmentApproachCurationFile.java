package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class RelevantTreatmentApproachCurationFile {

    private static final String FIELD_DELIMITER = "\t";

    @NotNull
    public static List<RelevantTreatmentApprochCurationEntry> read(@NotNull String drugClassCurationTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(drugClassCurationTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static List<RelevantTreatmentApprochCurationEntry> fromLines(@NotNull List<String> lines) {
        List<RelevantTreatmentApprochCurationEntry> curatedEntries = Lists.newArrayList();
        for (String line : lines) {
            curatedEntries.add(fromLine(line));
        }
        return curatedEntries;
    }

    @NotNull
    private static RelevantTreatmentApprochCurationEntry fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableRelevantTreatmentApprochCurationEntry.builder()
                .curationType(RelevantTreatmentApproachCurationType.valueOf(values[0]))
                .curationKey(ImmutableRelevantTreatmentApprochCurationEntryKey.builder()
                        .treatment(values[1])
                        .treatmentApproach(values[2])
                        .event(values[3])
                        .level(EvidenceLevel.valueOf(values[4]))
                        .direction(EvidenceDirection.valueOf(values[5]))
                        .build())
                .curatedtreatmentApproach(values.length == 5 ? values[6] : Strings.EMPTY)
                .build();
    }
}