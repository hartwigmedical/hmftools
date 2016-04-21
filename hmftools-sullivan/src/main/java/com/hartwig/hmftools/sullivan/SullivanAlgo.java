package com.hartwig.hmftools.sullivan;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.util.*;

public class SullivanAlgo {

    private static final Logger LOGGER = LogManager.getLogger(SullivanAlgo.class);

    @NotNull
    private final FileNameConverter originalFileNameConverter;

    public SullivanAlgo(@NotNull FileNameConverter originalFileNameConverter) {
        this.originalFileNameConverter = originalFileNameConverter;
    }

    public boolean runSullivanAlgo(@NotNull Iterable<String> fastqPaths,
                                   @NotNull String recreatedFastqPath, int numRecords) {
        boolean success = true;

        for (String fastqPath : fastqPaths) {
            success = success && runSullivanOnTwoDirectories(fastqPath, recreatedFastqPath, numRecords);
        }

        return success;
    }

    private boolean runSullivanOnTwoDirectories(@NotNull String originalFastqPath,
                                                @NotNull String recreatedFastqPath,
                                                int numRecords) {
        File originalPath = new File(originalFastqPath);
        File recreatedPath = new File(recreatedFastqPath);
        assert originalPath.isDirectory() && recreatedPath.isDirectory();

        File[] originalFiles = originalPath.listFiles();

        assert originalFiles != null;
        File[] recreatedFiles = new File[originalFiles.length];

        for (int i = 0; i < originalFiles.length; i++) {
            recreatedFiles[i] = new File(recreatedFastqPath + File.separator +
                    originalFileNameConverter.apply(originalFiles[i].getName()));
            LOGGER.info("Mapped " + originalFiles[i].getPath() + " to " + recreatedFiles[i].getPath());
        }

        boolean success = true;
        for (int i = 0; i < originalFiles.length; i++) {
            success = success && runSullivanOnTwoFiles(originalFiles[i], recreatedFiles[i], numRecords);
        }

        return success;
    }

    private static boolean runSullivanOnTwoFiles(@NotNull File originalFastqFile, @NotNull File recreatedFastqFile,
                                                 int numRecords) {
        assert originalFastqFile.isFile() && recreatedFastqFile.isFile();

        LOGGER.info("Start reading original fastq file from " + originalFastqFile.getPath());

        String refHeader = referenceHeader(originalFastqFile);
        if (refHeader == null) {
            LOGGER.warn("No ref header could be isolated from fastq file on " + originalFastqFile.getName());
            return false;
        } else {
            LOGGER.info("Generated ref header: " + refHeader);
        }

        Map<FastqHeaderKey, FastqRecord> originalFastq = mapOriginalFastqFile(originalFastqFile, refHeader, numRecords);
        int originalSize = originalFastq.size();
        LOGGER.info("Finished reading original fastq file. Created " + originalSize + " records.");

        FastqReader recreatedFastqReader = new FastqReader(recreatedFastqFile);
        FastqHeaderNormalizer recreatedNormalizer = new RecreatedFastqHeaderNormalizer();

        boolean success = true;
        int recordCount = 0;

        LOGGER.info("Start mapping process from recreated to original fastq");
        for (FastqRecord recreatedRecord : recreatedFastqReader) {
            FastqHeader header = FastqHeader.parseFromFastqRecord(recreatedRecord, recreatedNormalizer);
            if (!header.reference().equals(refHeader)) {
                LOGGER.warn("  Invalid header in recreated fastq file. Record = " + recreatedRecord);
            } else {
                FastqRecord originalMatch = originalFastq.get(header.key());
                if (originalMatch != null) {
                    if (!originalMatch.getReadString().equals(recreatedRecord.getReadString()) ||
                            !originalMatch.getBaseQualityString().equals(recreatedRecord.getBaseQualityString())) {
                        LOGGER.warn("  Mismatch between original and recreated fastq on record: " + recreatedRecord);
                        success = false;
                    }
                    originalFastq.remove(header.key());
                }
            }
            recordCount++;
            if (recordCount % 1E7 == 0) {
                int recordsFound = originalSize - originalFastq.size();
                LOGGER.info("  Finished mapping " + recordCount + " records. Found " + recordsFound + " original records");
            }
        }
        int recordsFound = originalSize - originalFastq.size();
        LOGGER.info("  Finished mapping " + recordCount + " records. Found " + recordsFound + " original records");

        LOGGER.info("Finished mapping records. " + originalFastq.size() + " unmapped records remaining in original fastq");

        return success && originalFastq.size() == 0;
    }

    @Nullable
    private static String referenceHeader(@NotNull File originalFastqFile) {
        FastqReader fastqReader = new FastqReader(originalFastqFile);
        if (fastqReader.hasNext()) {
            FastqHeader header = FastqHeader.parseFromFastqRecord(fastqReader.next(), new OriginalFastqHeaderNormalizer());
            fastqReader.close();
            return header.reference();
        }

        return null;
    }

    @NotNull
    private static Map<FastqHeaderKey, FastqRecord> mapOriginalFastqFile(
            @NotNull File file, @NotNull String refHeader, int numRecords) {
        FastqHeaderNormalizer normalizer = new OriginalFastqHeaderNormalizer();
        Map<FastqHeaderKey, FastqRecord> records = new HashMap<FastqHeaderKey, FastqRecord>(numRecords);

        FastqReader fastqReader = new FastqReader(file);

        while (fastqReader.hasNext() && records.size() < numRecords) {
            FastqRecord record = fastqReader.next();
            FastqHeader header = FastqHeader.parseFromFastqRecord(record, normalizer);

            if (records.containsKey(header.key())) {
                LOGGER.warn("  Duplicate record found: " + record);
            }

            records.put(header.key(), record);

            if (!header.reference().equals(refHeader)) {
                LOGGER.warn("  Header mismatch with reference header: " + header.reference());
            }

            if (records.size() % 1E7 == 0) {
                LOGGER.info("  Finished reading " + records.size() + " records");
            }
        }

        fastqReader.close();
        return records;
    }
}
