package com.hartwig.hmftools.common.amber.qc;

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

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".amber.qc";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
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
            String uniparentalDisomy = getValue(lines, UNIPARENTAL_DISOMY, null);

            if(uniparentalDisomy.equals(UNIPARENTAL_DISOMY_NONE))
                uniparentalDisomy = null;

            return ImmutableAmberQC.builder()
                    .contamination(Double.parseDouble(getValue(lines, CONTAMINATION, "0")))
                    .consanguinityProportion(Double.parseDouble(getValue(lines, CONSANGUINITY_PROPORTION, "0")))
                    .uniparentalDisomy(uniparentalDisomy)
                    .build();
        }
        catch(Exception e)
        {
            throw new IOException(String.format("Unable to parse amber qc file with %s lines.", lines.size()));
        }
    }

    private static String getValue(final List<String> lines, final String field, final String defaultValue)
    {
        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);
            if(values[0].equals(field))
                return values[1];
        }

        return defaultValue;
    }

    @NotNull
    private static List<String> toLines(@NotNull final AmberQC check)
    {
        final List<String> result = Lists.newArrayList();

        result.add(QC_STATUS + DELIMITER + check.status());
        result.add(CONTAMINATION + DELIMITER + FORMAT.format(check.contamination()));
        result.add(CONSANGUINITY_PROPORTION + DELIMITER + FORMAT.format(check.consanguinityProportion()));
        result.add(UNIPARENTAL_DISOMY + DELIMITER + (check.uniparentalDisomy() != null ? check.uniparentalDisomy() : UNIPARENTAL_DISOMY_NONE));

        return result;
    }
}
