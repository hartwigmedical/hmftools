package com.hartwig.hmftools.common.purple.purity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.qc.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;

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
    static PurpleQC fromLines(@NotNull final List<String> lines) {
        final ImmutablePurpleQC.Builder builder = ImmutablePurpleQC.builder().amberMeanDepth(0);

        if (lines.get(1).startsWith("SegmentPass")) {
            // PURPLE 2.47
            builder.unsupportedCopyNumberSegments(Integer.parseInt(getValue(lines.get(5))))
                    .amberGender(Gender.valueOf(getValue(lines.get(7))))
                    .cobaltGender(Gender.valueOf(getValue(lines.get(8))))
                    .deletedGenes(Integer.parseInt(getValue(lines.get(9))))
                    .copyNumberSegments(0)
                    .purity(1)
                    .method(FittedPurityMethod.NORMAL)
                    .contamination(0);
        } else {
            builder.method(FittedPurityMethod.valueOf(getValue(lines.get(1))))
                    .copyNumberSegments(Integer.parseInt(getValue(lines.get(2))))
                    .unsupportedCopyNumberSegments(Integer.parseInt(getValue(lines.get(3))))
                    .purity(Double.parseDouble(getValue(lines.get(4))))
                    .amberGender(Gender.valueOf(getValue(lines.get(5))))
                    .cobaltGender(Gender.valueOf(getValue(lines.get(6))))
                    .deletedGenes(Integer.parseInt(getValue(lines.get(7))))
                    .contamination(Double.parseDouble(getValue(lines.get(8))))
                    .germlineAberrations(GermlineAberration.fromString(getValue(lines.get(9))));
        }

        if (lines.size() > 10) {
            builder.amberMeanDepth(Integer.parseInt(getValue(lines.get(10))));
        }

        return builder.build();
    }

    @NotNull
    private static String getValue(@NotNull String line) {
        return line.split(DELIMITER)[1];
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final PurpleQC check) {
        final List<String> result = Lists.newArrayList();

        result.add("QCStatus" + DELIMITER + check.toString());
        result.add("Method" + DELIMITER + check.method());
        result.add("CopyNumberSegments" + DELIMITER + check.copyNumberSegments());
        result.add("UnsupportedCopyNumberSegments" + DELIMITER + check.unsupportedCopyNumberSegments());
        result.add("Purity" + DELIMITER + FORMATTER.format(check.purity()));
        result.add("AmberGender" + DELIMITER + check.amberGender());
        result.add("CobaltGender" + DELIMITER + check.cobaltGender());
        result.add("DeletedGenes" + DELIMITER + check.deletedGenes());
        result.add("Contamination" + DELIMITER + check.contamination());
        result.add("GermlineAberrations" + DELIMITER + GermlineAberration.toString(check.germlineAberrations()));
        result.add("AmberMeanDepth" + DELIMITER + check.amberMeanDepth());
        return result;
    }
}
