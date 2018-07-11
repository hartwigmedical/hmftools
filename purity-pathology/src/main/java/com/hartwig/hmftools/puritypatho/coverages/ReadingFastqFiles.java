package com.hartwig.hmftools.puritypatho.coverages;

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
      //  final FastqReader fastqReader = new FastqReader(fastqFile);

  //      for ( FastqRecord fastqRecord: fastqReader) {
    //        LOGGER.info(fastqRecord.getReadName());
      //  }

    }
}
