package com.hartwig.hmftools.puritypatho.coverages;

import java.io.File;
import java.io.IOException;

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
        for (String file:fileList) {
            if (file.endsWith(".fastq.gz")) {
                LOGGER.info(".fastq.gz file");
            } else if (file.endsWith(".fastq")) {
                LOGGER.info(".fastq file");
            }
            LOGGER.info("Processing file: " + file);
        }
    }
}
