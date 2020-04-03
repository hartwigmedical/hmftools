package com.hartwig.hmftools.sage.quality;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class QualityFile {

    private static final String DELIMITER = ",";

    public static void write(@NotNull final String filename, @NotNull final Collection<QualityCount> counts) throws IOException {
        Files.write(new File(filename).toPath(), toLines(counts));
    }

    @NotNull
    private static List<String> toLines(@NotNull final Collection<QualityCount> bafs) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        bafs.stream().map(QualityFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String toString(@NotNull final QualityCount baf) {
        return (char) baf.alt() + DELIMITER + (char) baf.ref() + DELIMITER + baf.qual() + DELIMITER + baf.firstOfPair() + DELIMITER
                + baf.readIndex() + DELIMITER + new String(baf.trinucleotideContext()) + DELIMITER + baf.count();

    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("alt")
                .add("ref")
                .add("qual")
                .add("firstOfPair")
                .add("readIndex")
                .add("trinucleotideContext")
                .add("count")
                .toString();
    }

}
