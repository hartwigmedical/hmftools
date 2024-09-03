package com.hartwig.hmftools.common.metrics;

import static com.hartwig.hmftools.common.metrics.BamMetricsSummary.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BamFlagStats
{
    public abstract long uniqueReadCount();
    public abstract long secondaryCount();
    public abstract long supplementaryCount();
    public abstract double duplicateProportion();
    public abstract double mappedProportion();
    public abstract double pairedInSequencingProportion();
    public abstract double properlyPairedProportion();
    public abstract double withItselfAndMateMappedProportion();
    public abstract double singletonProportion();

    public static final String FILE_EXTENSION = ".flag_counts.tsv";

    public static String generateFilename(final String basePath, final String sampleId)
    {
        return checkAddDirSeparator(basePath) + sampleId + BAM_METRICS_FILE_ID + FILE_EXTENSION;
    }

    public static final String FLAGSTAT_TOTAL = "in total (QC-passed reads + QC-failed reads)";
    public static final String FLAGSTAT_TOTAL_SHORT = "total";
    public static final String FLAGSTAT_PRIMARY = "primary";
    public static final String FLAGSTAT_SECONDARY = "secondary";
    public static final String FLAGSTAT_SUPPLEMENTARY = "supplementary";
    public static final String FLAGSTAT_DUPLICATE = "duplicates";
    public static final String FLAGSTAT_PRIMARY_DUPLICATE = "primary duplicates";
    public static final String FLAGSTAT_MAPPED = "mapped";
    public static final String FLAGSTAT_PRIMARY_MAPPED = "primary mapped";
    public static final String FLAGSTAT_PAIRED = "paired in sequencing";
    public static final String FLAGSTAT_READ1 = "read1";
    public static final String FLAGSTAT_READ2 = "read2";
    public static final String FLAGSTAT_PROPER_PAIR = "properly paired";
    public static final String FLAGSTAT_BOTH_MAPPED = "with itself and mate mapped";
    public static final String FLAGSTAT_SINGLE = "singletons";
    public static final String FLAGSTAT_MATE_DIFF_CHR = "with mate mapped to a different chr";
    public static final String FLAGSTAT_MATE_DIFF_CHR_MAPQ_5 = "with mate mapped to a different chr (mapQ>=5)";

    public static BamFlagStats read(final String flagstatPath) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(flagstatPath).toPath());
        String total = valueBySubstring(lines, FLAGSTAT_TOTAL_SHORT);
        String secondary = valueBySubstring(lines, FLAGSTAT_SECONDARY);
        String supplementary = valueBySubstring(lines, FLAGSTAT_SUPPLEMENTARY);
        String duplicates = valueBySubstring(lines, FLAGSTAT_DUPLICATE);
        String mapped = valueBySubstring(lines, "mapped (");
        String pairedInSequencing = valueBySubstring(lines, FLAGSTAT_PAIRED);
        String properlyPaired = valueBySubstring(lines, FLAGSTAT_PROPER_PAIR);
        String withItselfAndMateMapped = valueBySubstring(lines, FLAGSTAT_BOTH_MAPPED);
        String singletons = valueBySubstring(lines, FLAGSTAT_SINGLE);

        if(anyNull(
                total, secondary, supplementary, duplicates, mapped, pairedInSequencing, properlyPaired, withItselfAndMateMapped, singletons))
        {
            throw new IOException("Unable to parse flagstat file correctly");
        }

        long totalReadCount = Long.parseLong(total);
        long secondaryCount = Long.parseLong(secondary);
        long supplementaryCount = Long.parseLong(supplementary);

        // The total read count in flagstats double-counts secondary and supplementary reads, so need to remove to get unique read.
        long uniqueReadCount = totalReadCount - secondaryCount - supplementaryCount;

        // Paired-in-sequencing and properly-paired should be calculated over unique reads
        return ImmutableBamFlagStats.builder()
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
