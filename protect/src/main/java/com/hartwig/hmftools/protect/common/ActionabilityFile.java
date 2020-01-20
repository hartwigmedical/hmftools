package com.hartwig.hmftools.protect.common;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.protect.actionability.ActionabilitySource;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;
import com.hartwig.hmftools.protect.actionability.EvidenceLevel;
import com.hartwig.hmftools.protect.actionability.EvidenceScope;
import com.hartwig.hmftools.protect.actionability.ImmutableEvidenceItem;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActionabilityFile {

    private static final String DELIMITER = "\t";
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityFile.class);


    private ActionabilityFile() {

    }

    @NotNull
    public static List<EvidenceItem> read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    @NotNull
    private static List<EvidenceItem> fromLines(@NotNull final List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(ActionabilityFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static EvidenceItem fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        LOGGER.info(values[0]);
        final ImmutableEvidenceItem.Builder builder = ImmutableEvidenceItem.builder()
                .event(values[0])
                .source(ActionabilitySource.fromString(values[1]))
                .reference(values[2])
                .drug(values[3])
                .drugsType(values[4])
                .level(EvidenceLevel.fromString(values[5]))
                .response(values[6])
                .isOnLabel(Boolean.parseBoolean(values[7]))
                .cancerType(values[8])
                .scope((EvidenceScope.valueOf(values[9])));

        return builder.build();
    }
}
