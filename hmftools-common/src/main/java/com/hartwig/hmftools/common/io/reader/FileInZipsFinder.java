package com.hartwig.hmftools.common.io.reader;

import static java.util.stream.Collectors.toList;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.function.Predicate;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.jetbrains.annotations.NotNull;

@FunctionalInterface
interface FileInZipsFinder {

    String FILE_NOT_FOUND = "File %s not found in %s";

    @NotNull
    List<? extends ZipEntry> findFileInZip(@NotNull ZipFile zipFile, @NotNull String fileNameInZip)
            throws IOException;

    @NotNull
    static FileInZipsFinder build() {
        return (zipFile, fileNameInZip) -> {
            final Predicate<ZipEntry> isFile = zipEntry -> !zipEntry.isDirectory();
            final Predicate<ZipEntry> isFileName = zipEntry -> zipEntry.getName().contains(fileNameInZip);
            final List<? extends ZipEntry> fileInZip = zipFile.stream().filter(isFile.and(isFileName)).collect(
                    toList());
            if (fileInZip.isEmpty()) {
                throw new FileNotFoundException(String.format(FILE_NOT_FOUND, fileNameInZip, zipFile.getName()));
            }
            return fileInZip;
        };
    }
}
