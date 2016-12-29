package com.hartwig.healthchecker.common.io.reader;

import static java.util.stream.Collectors.toList;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import com.hartwig.healthchecker.common.exception.EmptyFileException;
import com.hartwig.healthchecker.common.exception.HealthChecksException;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
interface FileInZipsReader {

    @NotNull
    List<String> readLines(@NotNull String zipPath, @NotNull String fileNameInZip)
            throws IOException, HealthChecksException;

    @NotNull
    static FileInZipsReader build() {
        return (zipPath, fileNameInZip) -> {
            final List<String> fileLines = read(zipPath, fileNameInZip);
            if (fileLines.isEmpty()) {
                throw new EmptyFileException(fileNameInZip, zipPath);
            }
            return fileLines;
        };
    }

    @NotNull
    static List<String> read(@NotNull final String zipPath, @NotNull final String fileNameInZip) throws IOException {
        try (ZipFile zipFile = new ZipFile(zipPath)) {
            final List<? extends ZipEntry> fileEntryInZip = FileInZipsFinder.build().findFileInZip(zipFile,
                    fileNameInZip);
            return fileEntryInZip.stream().map(
                    zipElement -> ZipEntryReader.build().readZipElement(zipFile, zipElement).collect(
                            toList())).flatMap(Collection::stream).collect(toList());
        }
    }
}
