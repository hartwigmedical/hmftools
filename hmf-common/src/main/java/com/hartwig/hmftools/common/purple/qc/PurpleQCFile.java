package com.hartwig.hmftools.common.purple.qc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.io.exception.MalformedFileException;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFile {

    private static final DecimalFormat FORMATTER = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".purple.qc";

    private PurpleQCFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static PurpleQC read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final PurpleQC check) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    @NotNull
    @VisibleForTesting
    static PurpleQC fromLines(@NotNull final List<String> lines) throws MalformedFileException {
        try {
            return ImmutablePurpleQC.builder()
                    .unsupportedSegments(Integer.valueOf(getValue(lines.get(5))))
                    .ploidy(Double.valueOf(getValue(lines.get(6))))
                    .amberGender(Gender.valueOf(getValue(lines.get(7))))
                    .cobaltGender(Gender.valueOf(getValue(lines.get(8))))
                    .deletedGenes(Integer.valueOf(getValue(lines.get(9))))
                    .build();
        } catch (Exception e) {
            throw new MalformedFileException("Unable to parse purple qc file.");
        }
    }

    @NotNull
    private static String getValue(@NotNull String line) {
        return line.split(DELIMITER)[1];
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final PurpleQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.status());
        result.add("SegmentPass" + DELIMITER + check.segmentPass());
        result.add("GenderPass" + DELIMITER + check.genderPass());
        result.add("DeletedGenesPass" + DELIMITER + check.deletedGenesPass());
        result.add("SegmentScore" + DELIMITER + check.segmentScore());
        result.add("UnsupportedSegments" + DELIMITER + check.unsupportedSegments());
        result.add("Ploidy" + DELIMITER + FORMATTER.format(check.ploidy()));
        result.add("AmberGender" + DELIMITER + check.amberGender());
        result.add("CobaltGender" + DELIMITER + check.cobaltGender());
        result.add("DeletedGenes" + DELIMITER + check.deletedGenes());
        return result;
    }
}
