package com.hartwig.hmftools.healthchecker.io.reader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import com.hartwig.hmftools.healthchecker.exception.HealthChecksException;

import org.junit.Test;

import mockit.Mocked;
import mockit.NonStrictExpectations;

public class ZipEntryReaderTest {

    private static final String DATA = "Line 1";

    @Test
    public void readEntry(@Mocked final ZipFile zipFile, @Mocked final ZipEntry zipEntry,
                    @Mocked final InputStream inputStream, @Mocked final BufferedReader reader,
                    @Mocked final InputStreamReader inputStreamReader) throws IOException, HealthChecksException {
        new NonStrictExpectations() {
            {
                zipFile.getInputStream(zipEntry);
                result = inputStream;

                new InputStreamReader(inputStream);
                result = inputStreamReader;

                new BufferedReader(inputStreamReader);
                result = reader;

                reader.lines();
                returns(Stream.of(DATA, "", null));
            }
        };
        final Stream<String> zipStream = ZipEntryReader.build().readZipElement(zipFile, zipEntry);
        final List<String> lines = zipStream.collect(Collectors.toList());
        assertNotNull("zip entry is null", zipStream);
        assertEquals("Wrong num of lines", 1, lines.size());
    }

    @Test
    public void readEntryIoException(@Mocked final ZipFile zipFile, @Mocked final ZipEntry zipEntry,
                    @Mocked final InputStream inputStream, @Mocked final BufferedReader reader,
                    @Mocked final InputStreamReader inputStreamReader) throws IOException, HealthChecksException {
        new NonStrictExpectations() {
            {
                zipFile.getInputStream(zipEntry);
                result = new IOException();
            }
        };
        final Stream<String> zipStream = ZipEntryReader.build().readZipElement(zipFile, zipEntry);
        final List<String> lines = zipStream.collect(Collectors.toList());
        assertEquals("Wrong num of lines", 0, lines.size());
    }
}
