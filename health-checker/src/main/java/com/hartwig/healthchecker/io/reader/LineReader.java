package com.hartwig.healthchecker.io.reader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.healthchecker.exception.HealthChecksException;
import com.hartwig.healthchecker.exception.LineNotFoundException;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
public interface LineReader {

    @NotNull
    List<String> readLines(@NotNull Path filePath, @NotNull Predicate<String> filter)
            throws IOException, HealthChecksException;

    @NotNull
    static LineReader build() {
        return (filePath, filter) -> {
            final List<String> searchedLines = read(filePath, filter);
            if (searchedLines.isEmpty()) {
                throw new LineNotFoundException(filePath.toString(), filter.toString());
            }
            return searchedLines;
        };
    }

    @NotNull
    static List<String> read(@NotNull final Path filePath, @NotNull final Predicate<String> filter)
            throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(filePath.toString()))) {
            return lines.filter(filter).collect(Collectors.toList());
        }
    }
}
