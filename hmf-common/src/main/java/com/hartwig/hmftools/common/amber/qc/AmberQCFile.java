package com.hartwig.hmftools.common.amber.qc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.MalformedFileException;

import org.jetbrains.annotations.NotNull;

public class AmberQCFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".amber.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static AmberQC read(@NotNull final String filename) throws IOException, MalformedFileException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final AmberQC check) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    @NotNull
    static AmberQC fromLines(@NotNull final List<String> lines) throws IOException, MalformedFileException {
        try {
            return ImmutableAmberQC.builder()
                    .meanBAF(Double.valueOf(getValue(lines.get(2))))
                    .build();
        } catch (Exception e) {
            throw new MalformedFileException("Unable to parse amber qc file.");
        }
    }

    @NotNull
    private static String getValue(@NotNull String line) {
        return line.split(DELIMITER)[1];
    }

    @NotNull
    static List<String> toLines(@NotNull final AmberQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.status());
        result.add("MeanBAF" + DELIMITER + FORMAT.format(check.meanBAF()));

        return result;
    }
}
