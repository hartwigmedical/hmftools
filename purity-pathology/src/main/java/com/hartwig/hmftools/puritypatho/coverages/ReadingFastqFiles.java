package com.hartwig.hmftools.puritypatho.coverages;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReadingFastqFiles {
    private static final Logger LOGGER = LogManager.getLogger(ReadingFastqFiles.class);

    @NotNull
    public static void readingFiles(@NotNull String fastqDir) throws IOException {
        LOGGER.info("fastqDir: " + fastqDir);
        File folder = new File(fastqDir);
        String [] fileList = folder.list();
        final int size = 1048576;
        for (String file:fileList) {
            final InputStream inputStream;
            if (file.endsWith(".fastq.gz")) {
                inputStream = new GZIPInputStream(new FileInputStream(new File(fastqDir + file)), size);
                LOGGER.info("Processing file: " + fastqDir + file);
            } else if (file.endsWith(".fastq")) {
                inputStream = new FileInputStream(new File(fastqDir + file));
                LOGGER.info("Processing file: " + fastqDir + file);
            } else {
                LOGGER.info("Unknown file format");
            }
        }
    }
}
