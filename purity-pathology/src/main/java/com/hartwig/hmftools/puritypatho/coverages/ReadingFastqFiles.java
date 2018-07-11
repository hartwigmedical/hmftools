package com.hartwig.hmftools.puritypatho.coverages;

import java.io.File;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class ReadingFastqFiles {
    private static final Logger LOGGER = LogManager.getLogger(ReadingFastqFiles.class);

    @NotNull
    public static void readingFiles(@NotNull String fastqDir) {
        LOGGER.info("fastqDir: " + fastqDir);
        File folder = new File(fastqDir);
        String [] fileList = folder.list();
        for (String fileName:fileList) {
            LOGGER.info("fileName: " + fileName);

        }
    }
}
