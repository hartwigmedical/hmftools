package com.hartwig.hmftools.common.hotspot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class HotspotEvidenceFile {

    static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".hotspot.evidence";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull final List<HotspotEvidence> evidence) throws IOException {
        Files.write(new File(filename).toPath(), toLines(evidence));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<HotspotEvidence> evidence) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        evidence.stream().map(HotspotEvidenceFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("Chromosome")
                .add("Position")
                .add("Ref")
                .add("Alt")
                .add("Type")
                .add("TumorEvidence")
                .add("TumorReads")
                .add("NormalEvidence")
                .add("NormalReads")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final HotspotEvidence evidence) {
        return new StringJoiner(DELIMITER).add(evidence.chromosome())
                .add(String.valueOf(evidence.position()))
                .add(evidence.ref())
                .add(evidence.alt())
                .add(evidence.type().toString())
                .add(String.valueOf(evidence.tumorEvidence()))
                .add(String.valueOf(evidence.tumorReads()))
                .add(String.valueOf(evidence.normalEvidence()))
                .add(String.valueOf(evidence.normalReads()))
                .toString();
    }

}
