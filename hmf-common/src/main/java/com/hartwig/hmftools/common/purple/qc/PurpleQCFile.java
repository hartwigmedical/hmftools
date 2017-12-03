package com.hartwig.hmftools.common.purple.qc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.MalformedFileException;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public class PurpleQCFile {

    private static final NumberFormat formatter = new DecimalFormat("#0.00");
    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".purple.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static PurpleQC read(@NotNull final String filename) throws IOException, MalformedFileException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final PurpleQC check) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    @NotNull
    static PurpleQC fromLines(@NotNull final List<String> lines) throws IOException, MalformedFileException {
        try {
            return ImmutablePurpleQC.builder()
                    .unsupportedSegments(Integer.valueOf(getValue(lines.get(4))))
                    .ploidy(Double.valueOf(getValue(lines.get(5))))
                    .amberGender(Gender.valueOf(getValue(lines.get(6))))
                    .cobaltGender(Gender.valueOf(getValue(lines.get(7))))
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
    static List<String> toLines(@NotNull final PurpleQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.status());
        result.add("SegmentPass" + DELIMITER + check.segmentPass());
        result.add("GenderPass" + DELIMITER + check.genderPass());
        result.add("SegmentScore" + DELIMITER + check.segmentScore());
        result.add("UnsupportedSegments" + DELIMITER + check.unsupportedSegments());
        result.add("Ploidy" + DELIMITER + formatter.format(check.ploidy()));
        result.add("AmberGender" + DELIMITER + check.amberGender());
        result.add("CobaltGender" + DELIMITER + check.cobaltGender());

        return result;
    }
}
