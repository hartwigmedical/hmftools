package com.hartwig.hmftools.healthchecker.io.reader;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

@FunctionalInterface
interface ZipEntryReader {

    String ERROR_MSG = "Error occurred when reading file. Will return empty stream. Error -> %s";

    Logger LOGGER = LogManager.getLogger(ZipEntryReader.class);

    @NotNull
    Stream<String> readZipElement(@NotNull ZipFile zipFile, @NotNull ZipEntry zipEntry);

    @NotNull
    static ZipEntryReader build() {
        return (zipFile, zipEntry) -> {
            Stream<String> readLines = Stream.empty();
            try {
                final InputStream inputStream = zipFile.getInputStream(zipEntry);
                final InputStreamReader inputStreamReader = new InputStreamReader(inputStream);
                final BufferedReader reader = new BufferedReader(inputStreamReader);
                readLines = reader.lines();
            } catch (final IOException e) {
                LOGGER.error(String.format(ERROR_MSG, e.getMessage()));
            }
            return readLines.filter(line -> line != null && !line.isEmpty());
        };
    }
}
