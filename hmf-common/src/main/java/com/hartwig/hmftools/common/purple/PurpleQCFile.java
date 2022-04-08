package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.PurpleCommon.DELIMITER;
import static com.hartwig.hmftools.common.purple.PurpleQCStatus.fromString;
import static com.hartwig.hmftools.common.purple.purity.FittedPurityMethod.NORMAL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.getValue;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFile
{
    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String EXTENSION = ".purple.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static PurpleQC read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final PurpleQC check) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    private static final String QC_STATUS = "QCStatus";
    private static final String METHOD = "Method";
    private static final String CN_SEGMENTS = "CopyNumberSegments";
    private static final String UNSUPPORTED_CN_SEGMENTS = "UnsupportedCopyNumberSegments";
    private static final String PURITY = "Purity";
    private static final String AMBER_GENDER = "AmberGender";
    private static final String COBALT_GENDER = "CobaltGender";
    private static final String DELETED_GENES = "DeletedGenes";
    private static final String CONTAMINATION = "Contamination";
    private static final String GERMLINE_ABERRATIONS = "GermlineAberrations";
    private static final String AMBER_MEAN_DEPTH = "AmberMeanDepth";

    private static PurpleQC fromLines(@NotNull final List<String> lines)
    {
        final ImmutablePurpleQC.Builder builder = ImmutablePurpleQC.builder();

        String qcStatusValues = getValue(lines, QC_STATUS, "", DELIMITER);
        Set<PurpleQCStatus> statusSet = PurpleQCStatus.fromString(qcStatusValues);

        builder.method(FittedPurityMethod.valueOf(getValue(lines, METHOD, NORMAL.toString(), DELIMITER)))
                .status(statusSet)
                .copyNumberSegments(Integer.parseInt(getValue(lines, CN_SEGMENTS, "0", DELIMITER)))
                .unsupportedCopyNumberSegments(Integer.parseInt(getValue(lines, UNSUPPORTED_CN_SEGMENTS, "0", DELIMITER)))
                .purity(Double.parseDouble(getValue(lines, PURITY, "0", DELIMITER)))
                .amberGender(Gender.valueOf(getValue(lines, AMBER_GENDER, null, DELIMITER)))
                .cobaltGender(Gender.valueOf(getValue(lines, COBALT_GENDER, null, DELIMITER)))
                .deletedGenes(Integer.parseInt(getValue(lines, DELETED_GENES, "0", DELIMITER)))
                .contamination(Double.parseDouble(getValue(lines, CONTAMINATION, "0", DELIMITER)))
                .germlineAberrations(GermlineAberration.fromString(getValue(
                        lines, GERMLINE_ABERRATIONS, GermlineAberration.NONE.toString(), DELIMITER)))
                .amberMeanDepth(Integer.parseInt(getValue(lines, AMBER_MEAN_DEPTH, "0", DELIMITER)));

        return builder.build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final PurpleQC check)
    {
        final List<String> result = Lists.newArrayList();

        result.add(QC_STATUS + DELIMITER + check.toString());
        result.add(METHOD + DELIMITER + check.method());
        result.add(CN_SEGMENTS + DELIMITER + check.copyNumberSegments());
        result.add(UNSUPPORTED_CN_SEGMENTS + DELIMITER + check.unsupportedCopyNumberSegments());
        result.add(PURITY + DELIMITER + FORMAT.format(check.purity()));
        result.add(AMBER_GENDER + DELIMITER + check.amberGender());
        result.add(COBALT_GENDER + DELIMITER + check.cobaltGender());
        result.add(DELETED_GENES + DELIMITER + check.deletedGenes());
        result.add(CONTAMINATION + DELIMITER + check.contamination());
        result.add(GERMLINE_ABERRATIONS + DELIMITER + GermlineAberration.toString(check.germlineAberrations()));
        result.add(AMBER_MEAN_DEPTH + DELIMITER + check.amberMeanDepth());
        return result;
    }
}
