package com.hartwig.hmftools.common.lims;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Optional;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.utils.io.exception.EmptyFileException;
import com.hartwig.hmftools.common.utils.io.exception.MalformedFileException;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class LimsWideFile {

    private static final String DELIMITER = "\t";

    private LimsWideFile() {
    }


    @NotNull
    public static LimsWide read(@NotNull String filePath) throws IOException {
        ImmutableLimsWide.Builder builder = ImmutableLimsWide.builder();
        String[] values = findValuesLine(filePath).split(DELIMITER);

        builder.studyName(values[0]);
        builder.reportReceiverName(values[1]);
        builder.reportReceiverEmail(values[2]);

        return builder.build();
    }

    @NotNull
    public static LimsWide empty() {
        return ImmutableLimsWide.builder().studyName(Strings.EMPTY).reportReceiverName(Strings.EMPTY).reportReceiverEmail(Strings.EMPTY).build();
    }

    @NotNull
    private static String findValuesLine(@NotNull String filename) throws IOException {
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if (lines.isEmpty()) {
            throw new EmptyFileException(filename);
        }
        final int index = findHeaderLineIndex(lines);
        if (index >= lines.size()) {
            throw new MalformedFileException(String.format("No value line found after header line in contact file %s.", filename));
        }
        return lines.get(index + 1);
    }

    private static int findHeaderLineIndex(@NotNull final List<String> lines) throws MalformedFileException {
        final Optional<Integer> lineNumbers =
                IntStream.range(0, lines.size()).filter(index -> lines.get(index).contains("StudyName")).boxed().findFirst();
        if (!lineNumbers.isPresent()) {
            throw new MalformedFileException(String.format("Could not find header line in contact file with %s lines.", lines.size()));
        }
        return lineNumbers.get();
    }

}
