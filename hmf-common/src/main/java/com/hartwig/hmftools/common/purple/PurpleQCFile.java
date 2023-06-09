package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.FittedPurityMethod.NORMAL;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.getValue;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFile
{
    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String EXTENSION = ".purple.qc";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static PurpleQC read(final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(final String fileName, final PurpleQC check) throws IOException
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

    private static PurpleQC fromLines(final List<String> lines)
    {
        final ImmutablePurpleQC.Builder builder = ImmutablePurpleQC.builder();

        String qcStatusValues = getValue(lines, QC_STATUS, "", TSV_DELIM);
        Set<PurpleQCStatus> statusSet = PurpleQCStatus.fromString(qcStatusValues);

        builder.method(FittedPurityMethod.valueOf(getValue(lines, METHOD, NORMAL.toString(), TSV_DELIM)))
                .status(statusSet)
                .copyNumberSegments(Integer.parseInt(getValue(lines, CN_SEGMENTS, "0", TSV_DELIM)))
                .unsupportedCopyNumberSegments(Integer.parseInt(getValue(lines, UNSUPPORTED_CN_SEGMENTS, "0", TSV_DELIM)))
                .purity(Double.parseDouble(getValue(lines, PURITY, "0", TSV_DELIM)))
                .amberGender(Gender.valueOf(getValue(lines, AMBER_GENDER, null, TSV_DELIM)))
                .cobaltGender(Gender.valueOf(getValue(lines, COBALT_GENDER, null, TSV_DELIM)))
                .deletedGenes(Integer.parseInt(getValue(lines, DELETED_GENES, "0", TSV_DELIM)))
                .contamination(Double.parseDouble(getValue(lines, CONTAMINATION, "0", TSV_DELIM)))
                .germlineAberrations(GermlineAberration.fromString(getValue(
                        lines, GERMLINE_ABERRATIONS, GermlineAberration.NONE.toString(), TSV_DELIM)))
                .amberMeanDepth(Integer.parseInt(getValue(lines, AMBER_MEAN_DEPTH, "0", TSV_DELIM)));

        return builder.build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(final PurpleQC check)
    {
        final List<String> result = Lists.newArrayList();

        result.add(QC_STATUS + TSV_DELIM + check.toString());
        result.add(METHOD + TSV_DELIM + check.method());
        result.add(CN_SEGMENTS + TSV_DELIM + check.copyNumberSegments());
        result.add(UNSUPPORTED_CN_SEGMENTS + TSV_DELIM + check.unsupportedCopyNumberSegments());
        result.add(PURITY + TSV_DELIM + FORMAT.format(check.purity()));
        result.add(AMBER_GENDER + TSV_DELIM + check.amberGender());
        result.add(COBALT_GENDER + TSV_DELIM + check.cobaltGender());
        result.add(DELETED_GENES + TSV_DELIM + check.deletedGenes());
        result.add(CONTAMINATION + TSV_DELIM + check.contamination());
        result.add(GERMLINE_ABERRATIONS + TSV_DELIM + GermlineAberration.toString(check.germlineAberrations()));
        result.add(AMBER_MEAN_DEPTH + TSV_DELIM + check.amberMeanDepth());
        return result;
    }
}
