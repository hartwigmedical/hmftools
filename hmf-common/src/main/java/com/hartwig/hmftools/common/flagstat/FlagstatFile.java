package com.hartwig.hmftools.common.flagstat;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FlagstatFile
{
    public static final String FILE_EXTENSION = ".flagstat";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + FILE_EXTENSION;
    }

    public static Flagstat read(final String flagstatPath) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(flagstatPath).toPath());
        String total = valueBySubstring(lines, "total");
        String secondary = valueBySubstring(lines, "secondary");
        String supplementary = valueBySubstring(lines, "supplementary");
        String duplicates = valueBySubstring(lines, "duplicates");
        String mapped = valueBySubstring(lines, "mapped (");
        String pairedInSequencing = valueBySubstring(lines, "paired in sequencing");
        String properlyPaired = valueBySubstring(lines, "properly paired");
        String withItselfAndMateMapped = valueBySubstring(lines, "with itself and mate mapped");
        String singletons = valueBySubstring(lines, "singletons");
        if(anyNull(total,
                secondary,
                supplementary,
                duplicates,
                mapped,
                pairedInSequencing,
                properlyPaired,
                withItselfAndMateMapped,
                singletons))
        {
            throw new IOException("Unable to parse flagstat file correctly");
        }

        long totalReadCount = Long.parseLong(total);
        long secondaryCount = Long.parseLong(secondary);
        long supplementaryCount = Long.parseLong(supplementary);

        // The total read count in flagstats double-counts secondary and supplementary reads, so need to remove to get unique read.
        long uniqueReadCount = totalReadCount - secondaryCount - supplementaryCount;

        // Paired-in-sequencing and properly-paired should be calculated over unique reads.
        return ImmutableFlagstat.builder()
                .uniqueReadCount(uniqueReadCount)
                .secondaryCount(secondaryCount)
                .supplementaryCount(supplementaryCount)
                .duplicateProportion(divideTwoStrings(duplicates, total))
                .mappedProportion(divideTwoStrings(mapped, total))
                .pairedInSequencingProportion(divideTwoStrings(pairedInSequencing, String.valueOf(uniqueReadCount)))
                .properlyPairedProportion(divideTwoStrings(properlyPaired, String.valueOf(uniqueReadCount)))
                .withItselfAndMateMappedProportion(divideTwoStrings(withItselfAndMateMapped, total))
                .singletonProportion(divideTwoStrings(singletons, total))
                .build();
    }

    private static double divideTwoStrings(final String string1, final String string2)
    {
        long value1 = Long.parseLong(string1);
        long value2 = Long.parseLong(string2);
        return roundToSixDecimals((double)value1 / value2);
    }

    private static double roundToSixDecimals(double value)
    {
        BigDecimal bd = BigDecimal.valueOf(value);
        bd = bd.setScale(6, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    @Nullable
    private static String valueBySubstring(final List<String> lines, final String subString)
    {
        List<String> matchLines = Lists.newArrayList();
        for(String line : lines)
        {
            if(line.contains(subString))
            {
                matchLines.add(line);
            }
        }
        if(matchLines.size() >= 1)
        {
            return matchLines.get(0).split(" ")[0];
        }
        return null;
    }

    private static boolean anyNull(final Object... arguments)
    {
        for(Object object : arguments)
        {
            if(object == null)
            {
                return true;
            }
        }
        return false;
    }
}
