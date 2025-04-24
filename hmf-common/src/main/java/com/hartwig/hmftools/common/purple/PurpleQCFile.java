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

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

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
    private static final String LOH_PERCENT = "LohPercent";
    private static final String TINC_LEVEL = "TincLevel";
    private static final String CHIMERISM_PRESENT = "ChimerismPresent";
    private static final String CHIMERISM_PERCENTAGE = "ChimerismPercentage";

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
                .amberMeanDepth(Integer.parseInt(getValue(lines, AMBER_MEAN_DEPTH, "0", TSV_DELIM)))
                .lohPercent(Double.parseDouble(getValue(lines, LOH_PERCENT, "-1", TSV_DELIM)))
                .tincLevel(Double.parseDouble(getValue(lines, TINC_LEVEL, "0", TSV_DELIM)))
                .chimerismPresent(Boolean.parseBoolean(getValue(lines, CHIMERISM_PRESENT, "false", TSV_DELIM)))
                .chimerismPercentage(Double.parseDouble(getValue(lines, CHIMERISM_PERCENTAGE, "-1", TSV_DELIM)));

        return builder.build();
    }

    @VisibleForTesting
    static List<String> toLines(final PurpleQC purpleQC)
    {
        final List<String> result = Lists.newArrayList();

        result.add(QC_STATUS + TSV_DELIM + purpleQC.toString());
        result.add(METHOD + TSV_DELIM + purpleQC.method());
        result.add(CN_SEGMENTS + TSV_DELIM + purpleQC.copyNumberSegments());
        result.add(UNSUPPORTED_CN_SEGMENTS + TSV_DELIM + purpleQC.unsupportedCopyNumberSegments());
        result.add(PURITY + TSV_DELIM + FORMAT.format(purpleQC.purity()));
        result.add(AMBER_GENDER + TSV_DELIM + purpleQC.amberGender());
        result.add(COBALT_GENDER + TSV_DELIM + purpleQC.cobaltGender());
        result.add(DELETED_GENES + TSV_DELIM + purpleQC.deletedGenes());
        result.add(CONTAMINATION + TSV_DELIM + purpleQC.contamination());
        result.add(GERMLINE_ABERRATIONS + TSV_DELIM + GermlineAberration.toString(purpleQC.germlineAberrations()));
        result.add(AMBER_MEAN_DEPTH + TSV_DELIM + purpleQC.amberMeanDepth());
        result.add(LOH_PERCENT + TSV_DELIM + FORMAT.format(purpleQC.lohPercent()));
        result.add(TINC_LEVEL + TSV_DELIM + FORMAT.format(purpleQC.tincLevel()));
        result.add(CHIMERISM_PRESENT + TSV_DELIM + purpleQC.chimerismPresent());
        if(purpleQC.chimerismPercentage() != null)
        {
            result.add(CHIMERISM_PERCENTAGE + TSV_DELIM + FORMAT.format(purpleQC.chimerismPercentage()));
        }
        return result;
    }
}
