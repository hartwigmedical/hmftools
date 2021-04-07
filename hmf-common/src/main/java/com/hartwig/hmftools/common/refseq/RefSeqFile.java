package com.hartwig.hmftools.common.refseq;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class RefSeqFile {

    private static final String FIELD_DELIMITER = "\t";

    private RefSeqFile() {
    }

    @NotNull
    public static List<RefSeq> readingRefSeq(@NotNull String refSeqFile) throws IOException {
        List<RefSeq> refSeq = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(refSeqFile).toPath());
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            refSeq.add(fromLine(line));
        }

        return refSeq;
    }

    @NotNull
    private static RefSeq fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableRefSeq.builder().geneId(values[0]).transcriptId(values[1]).displayLabel(values[2]).dbPrimaryAcc(values[3]).build();
    }
}
