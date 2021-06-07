package com.hartwig.hmftools.sage.quality;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class QualityRecalibrationFile
{

    private static final String DELIMITER = "\t";
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    public static void write(@NotNull final String filename, @NotNull final Collection<QualityRecalibrationRecord> counts)
            throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(counts));
    }

    @NotNull
    private static List<String> toLines(@NotNull final Collection<QualityRecalibrationRecord> bafs)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        bafs.stream().map(QualityRecalibrationFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String toString(@NotNull final QualityRecalibrationRecord baf)
    {
        return (char) baf.key().alt() + DELIMITER + (char) baf.key().ref() + DELIMITER + new String(baf.key().trinucleotideContext())
                + DELIMITER + baf.count() + DELIMITER + baf.key().qual() + DELIMITER + FORMAT.format(baf.recalibratedQual());
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "").add("alt")
                .add("ref")
                .add("trinucleotideContext")
                .add("count")
                .add("originalQual")
                .add("recalibratedQual")
                .toString();
    }
}
