package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class QualityRecalibrationFile
{
    public static void write(final String filename, final List<QualityRecalibrationRecord> counts) throws IOException
    {
        Collections.sort(counts);
        Files.write(new File(filename).toPath(), toLines(counts));
    }

    private static List<String> toLines(final Collection<QualityRecalibrationRecord> bafs)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        bafs.stream().map(QualityRecalibrationFile::toString).forEach(lines::add);
        return lines;
    }

    private static String toString(final QualityRecalibrationRecord baf)
    {
        StringJoiner sj = new StringJoiner(DELIM);
        sj.add(String.valueOf((char)baf.Key.Alt));
        sj.add(String.valueOf((char)baf.Key.Ref));
        sj.add(new String(baf.Key.TrinucleotideContext));
        sj.add(String.valueOf(baf.Count));
        sj.add(String.valueOf(baf.Key.Quality));
        sj.add(String.format("%.2f", baf.RecalibratedQuality));
        return sj.toString();
    }

    private static String header()
    {
        return new StringJoiner(DELIM, "", "")
                .add("alt")
                .add("ref")
                .add("trinucleotideContext")
                .add("count")
                .add("originalQual")
                .add("recalibratedQual")
                .toString();
    }
}
