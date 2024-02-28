package com.hartwig.hmftools.common.basequal.jitter;

import java.io.File;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.StringUtils;
import org.jetbrains.annotations.NotNull;

//
public class MicrosatelliteSiteFile
{
    private static final String CHROMOSOME = "chromosome";
    private static final String START = "start";
    private static final String END = "end";
    private static final String UNIT = "unit";
    private static final String NUM_READS = "numReads";
    private static final String NUM_READS_REJECTED = "numReadsRejected";

    private static final String REAL_VARIANT = "realVariant";

    private static final String COUNT_p10 = "count+10";
    private static final String COUNT_p9 = "count+9";
    private static final String COUNT_p8 = "count+8";
    private static final String COUNT_p7 = "count+7";
    private static final String COUNT_p6 = "count+6";
    private static final String COUNT_p5 = "count+5";
    private static final String COUNT_p4 = "count+4";
    private static final String COUNT_p3 = "count+3";
    private static final String COUNT_p2 = "count+2";
    private static final String COUNT_p1 = "count+1";
    private static final String COUNT_p0 = "count+0";

    private static final String COUNT_m10 = "count-10";
    private static final String COUNT_m9 = "count-9";
    private static final String COUNT_m8 = "count-8";
    private static final String COUNT_m7 = "count-7";
    private static final String COUNT_m6 = "count-6";
    private static final String COUNT_m5 = "count-5";
    private static final String COUNT_m4 = "count-4";
    private static final String COUNT_m3 = "count-3";
    private static final String COUNT_m2 = "count-2";
    private static final String COUNT_m1 = "count-1";

    private static final String READ_REPEAT_LENGTHS = "readRepeatLengths";

    private static final String FILE_EXTENSION = ".repeat.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers)
    {
        Comparator<MicrosatelliteSiteAnalyser> comparator =
                Comparator.comparing((MicrosatelliteSiteAnalyser o) -> o.refGenomeMicrosatellite.chromosome())
                        .thenComparingInt(o -> o.refGenomeMicrosatellite.referenceStart())
                        .thenComparing(o -> o.refGenomeMicrosatellite.referenceEnd());

        // sort the bins
        List<MicrosatelliteSiteAnalyser> sortedAnalysers = microsatelliteSiteAnalysers.stream().sorted(comparator).collect(Collectors.toList());

        List<String> columns = List.of(CHROMOSOME, START, END, UNIT, NUM_READS, NUM_READS_REJECTED, REAL_VARIANT,
                COUNT_m10, COUNT_m9, COUNT_m8, COUNT_m7, COUNT_m6, COUNT_m5, COUNT_m4, COUNT_m3, COUNT_m2, COUNT_m1,
                COUNT_p0, COUNT_p1, COUNT_p2, COUNT_p3, COUNT_p4, COUNT_p5, COUNT_p6, COUNT_p7, COUNT_p8, COUNT_p9, COUNT_p10,
                READ_REPEAT_LENGTHS);

        // add the count columns

        DelimFileWriter.write(filename, columns, sortedAnalysers, (repeatAnalyser, row) ->
        {
            row.set(CHROMOSOME, repeatAnalyser.refGenomeMicrosatellite.chromosome());
            row.set(START, repeatAnalyser.refGenomeMicrosatellite.referenceStart());
            row.set(END, repeatAnalyser.refGenomeMicrosatellite.referenceEnd());
            row.set(UNIT,  repeatAnalyser.refGenomeMicrosatellite.unitString());
            row.set(NUM_READS, repeatAnalyser.getReadRepeatMatches().size());
            row.set(NUM_READS_REJECTED, repeatAnalyser.getReadRepeatMatches().stream().filter(o -> o.shouldDropRead).count());
            row.set(REAL_VARIANT, repeatAnalyser.isRealVariant(JitterAnalyserConstants.ALT_COUNT_FRACTION_INIT, JitterAnalyserConstants.ALT_COUNT_FRACTION_STEP,
                    JitterAnalyserConstants.MAX_REJECTED_READ_FRACTION));
            int refNumRepeat = repeatAnalyser.refGenomeMicrosatellite.numRepeat;
            row.set(COUNT_p0, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat));

            row.set(COUNT_p10, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 10));
            row.set(COUNT_p9, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 9));
            row.set(COUNT_p8, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 8));
            row.set(COUNT_p7, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 7));
            row.set(COUNT_p6, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 6));
            row.set(COUNT_p5, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 5));
            row.set(COUNT_p4, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 4));
            row.set(COUNT_p3, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 3));
            row.set(COUNT_p2, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 2));
            row.set(COUNT_p1, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat + 1));

            row.set(COUNT_m10, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 10));
            row.set(COUNT_m9, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 9));
            row.set(COUNT_m8, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 8));
            row.set(COUNT_m7, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 7));
            row.set(COUNT_m6, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 6));
            row.set(COUNT_m5, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 5));
            row.set(COUNT_m4, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 4));
            row.set(COUNT_m3, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 3));
            row.set(COUNT_m2, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 2));
            row.set(COUNT_m1, repeatAnalyser.getCountWithRepeatUnits(refNumRepeat - 1));

            row.set(READ_REPEAT_LENGTHS, getRepeatString(repeatAnalyser.getPassingReadRepeatMatches()));
        });
    }

    public static String getRepeatString(final List<MicrosatelliteRead> microsatelliteReads)
    {
        List<Integer> repeatLengths = microsatelliteReads.stream().filter(o -> !o.shouldDropRead).map(MicrosatelliteRead::readRepeatLength).collect(
                Collectors.toList());
        return StringUtils.join(repeatLengths, ",");
    }
}
