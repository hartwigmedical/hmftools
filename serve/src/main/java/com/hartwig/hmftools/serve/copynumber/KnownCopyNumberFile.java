package com.hartwig.hmftools.serve.copynumber;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class KnownCopyNumberFile {

    private static final String DELIMITER = "\t";

    private KnownCopyNumberFile() {
    }

    public static void write(@NotNull String copyNumberTsv, @NotNull List<KnownCopyNumber> copyNumbers) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(copyNumbers));
        Files.write(new File(copyNumberTsv).toPath(), lines);
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene").add("type").add("sources").toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull List<KnownCopyNumber> copyNumbers) {
        List<String> lines = Lists.newArrayList();
        for (KnownCopyNumber copyNumber : copyNumbers) {
            lines.add(toLine(copyNumber));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull KnownCopyNumber copyNumber) {
        return new StringJoiner(DELIMITER).add(copyNumber.gene())
                .add(copyNumber.type().toString().toLowerCase())
                .add(Knowledgebase.commaSeparatedSourceString(copyNumber.sources()))
                .toString();
    }
}
