package com.hartwig.hmftools.common.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class AmberQCFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String EXTENSION = ".amber.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    public static AmberQC read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull final AmberQC check) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(check));
    }

    private static final String QC_STATUS = "QCStatus";
    private static final String CONTAMINATION = "Contamination";
    private static final String CONSANGUINITY_PROPORTION = "ConsanguinityProportion";
    private static final String UNIPARENTAL_DISOMY = "UniparentalDisomy";
    private static final String UNIPARENTAL_DISOMY_NONE = "NONE";

    @NotNull
    private static AmberQC fromLines(final List<String> lines) throws IOException
    {
        try
        {
            String uniparentalDisomy = getValue(lines, UNIPARENTAL_DISOMY, null, TSV_DELIM);

            if(uniparentalDisomy != null && uniparentalDisomy.equals(UNIPARENTAL_DISOMY_NONE))
                uniparentalDisomy = null;

            return ImmutableAmberQC.builder()
                    .contamination(Double.parseDouble(getValue(lines, CONTAMINATION, "0", TSV_DELIM)))
                    .consanguinityProportion(Double.parseDouble(getValue(lines, CONSANGUINITY_PROPORTION, "0", TSV_DELIM)))
                    .uniparentalDisomy(uniparentalDisomy)
                    .build();
        }
        catch(Exception e)
        {
            throw new IOException(String.format("Unable to parse amber qc file with %s lines.", lines.size()));
        }
    }

    @NotNull
    private static List<String> toLines(@NotNull final AmberQC check)
    {
        final List<String> result = Lists.newArrayList();

        result.add(QC_STATUS + TSV_DELIM + check.status());
        result.add(CONTAMINATION + TSV_DELIM + FORMAT.format(check.contamination()));
        result.add(CONSANGUINITY_PROPORTION + TSV_DELIM + FORMAT.format(check.consanguinityProportion()));
        result.add(UNIPARENTAL_DISOMY + TSV_DELIM + (check.uniparentalDisomy() != null ? check.uniparentalDisomy() : UNIPARENTAL_DISOMY_NONE));

        return result;
    }
}
