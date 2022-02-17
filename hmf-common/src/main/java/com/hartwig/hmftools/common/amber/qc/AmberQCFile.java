package com.hartwig.hmftools.common.amber.qc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class AmberQCFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".amber.qc";

    private AmberQCFile() {
    }

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static AmberQC read(@NotNull final String filename) throws IOException {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final AmberQC check) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    @NotNull
    private static AmberQC fromLines(@NotNull final List<String> lines) throws IOException {
        try {
            return ImmutableAmberQC.builder()
                    .meanBAF(Double.parseDouble(getValue(lines.get(1))))
                    .contamination(lines.size() > 2 ? Double.parseDouble(getValue(lines.get(2))) : 0)
                    .consanguinityProportion(lines.size() > 3 ? Double.parseDouble(getValue(lines.get(3))) : 0)
                    .uniparentalDisomy(lines.size() > 4 ? lines.get(4) : null)
                    .build();
        } catch (Exception e) {
            throw new IOException(String.format("Unable to parse amber qc file with %s lines.", lines.size()));
        }
    }

    @NotNull
    private static String getValue(@NotNull String line) {
        return line.split(DELIMITER)[1];
    }

    @NotNull
    private static List<String> toLines(@NotNull final AmberQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.status());
        result.add("MeanBAF" + DELIMITER + FORMAT.format(check.meanBAF()));
        result.add("Contamination" + DELIMITER + FORMAT.format(check.contamination()));
        result.add("ConsanguinityProportion" + DELIMITER + FORMAT.format(check.consanguinityProportion()));
        result.add("UniparentalDisomy" + DELIMITER + (check.uniparentalDisomy() != null ? check.uniparentalDisomy() : "null"));

        return result;
    }
}
