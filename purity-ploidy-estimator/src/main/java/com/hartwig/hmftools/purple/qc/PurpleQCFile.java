package com.hartwig.hmftools.purple.qc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class PurpleQCFile {

    private static final NumberFormat formatter = new DecimalFormat("#0.00");
    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".purple.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String fileName, @NotNull final PurpleQC check) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    private static List<String> toLines(@NotNull final PurpleQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.status());
        result.add("SegmentPass" + DELIMITER + check.segmentPass());
        result.add("GenderPass" + DELIMITER + check.genderPass());
        result.add("SegmentScore" + DELIMITER + check.segmentScore());
        result.add("TrailingSegments" + DELIMITER + check.trailingSegments());
        result.add("RatioOnlySegments" + DELIMITER + check.ratioSegments());
        result.add("Ploidy" + DELIMITER + formatter.format(check.ploidy()));
        result.add("PurpleGender" + DELIMITER + check.purpleGender());
        result.add("CobaltGender" + DELIMITER + check.cobaltGender());

        return result;
    }

}
