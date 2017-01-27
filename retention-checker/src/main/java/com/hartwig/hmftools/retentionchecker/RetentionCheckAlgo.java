package com.hartwig.hmftools.retentionchecker;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

class RetentionCheckAlgo {

    private static final Logger LOGGER = LogManager.getLogger(RetentionCheckAlgo.class);

    @NotNull
    private final FileNameConverter originalFileNameConverter;

    RetentionCheckAlgo(@NotNull final FileNameConverter originalFileNameConverter) {
        this.originalFileNameConverter = originalFileNameConverter;
    }

    boolean runAlgo(@NotNull final Iterable<String> fastqPaths, @NotNull final String recreatedFastqPath,
            final int numRecords) {
        boolean success = true;

        for (final String fastqPath : fastqPaths) {
            success = success && runRetentionCheckOnTwoDirectories(fastqPath, recreatedFastqPath, numRecords);
        }

        return success;
    }

    private boolean runRetentionCheckOnTwoDirectories(@NotNull String originalFastqPath,
            @NotNull String recreatedFastqPath, int numRecords) {
        final File originalPath = new File(originalFastqPath);
        final File recreatedPath = new File(recreatedFastqPath);
        if (!originalPath.exists()) {
            LOGGER.warn("Original fastq path does not exist: " + originalFastqPath);
            return false;
        } else if (!recreatedPath.exists()) {
            LOGGER.warn("Recreated fastq path does not exist: " + recreatedFastqPath);
            return false;
        }

        assert originalPath.isDirectory() && recreatedPath.isDirectory();

        final File[] originalFiles = originalPath.listFiles((dir, name) -> name.indexOf(".fastq") > 0);

        assert originalFiles != null;
        final File[] recreatedFiles = new File[originalFiles.length];

        for (int i = 0; i < originalFiles.length; i++) {
            recreatedFiles[i] = new File(
                    recreatedFastqPath + File.separator + originalFileNameConverter.apply(originalFiles[i].getName()));
            LOGGER.info("Mapped " + originalFiles[i].getPath() + " to " + recreatedFiles[i].getPath());
        }

        boolean success = true;
        for (int i = 0; i < originalFiles.length; i++) {
            success = success && runRetentionCheckOnTwoFiles(originalFiles[i], recreatedFiles[i], numRecords);
        }

        return success;
    }

    private static boolean runRetentionCheckOnTwoFiles(@NotNull final File originalFastqFile,
            @NotNull final File recreatedFastqFile, final int numRecords) {
        assert originalFastqFile.isFile() && recreatedFastqFile.isFile();

        LOGGER.info("Start reading recreated fastq file from " + recreatedFastqFile.getPath());

        final String refHeader = referenceHeader(recreatedFastqFile);
        if (refHeader == null) {
            LOGGER.warn("No ref header could be isolated from fastq file on " + recreatedFastqFile.getName());
            return false;
        } else {
            LOGGER.info("Generated ref header: " + refHeader);
        }

        final Map<FastqHeaderKey, FastqRecord> recreatedFastq = mapRecreatedFastqFile(recreatedFastqFile, refHeader,
                numRecords);
        final int recreatedSize = recreatedFastq.size();
        LOGGER.info("Finished reading recreated fastq file. Created " + recreatedSize + " records.");

        final FastqReader originalFastqReader = new FastqReader(originalFastqFile);
        final FastqHeaderNormalizer originalNormalizer = new OriginalFastqHeaderNormalizer();

        boolean success = true;
        int recordCount = 0;

        LOGGER.info("Start mapping process from original to recreated fastq");
        for (final FastqRecord originalRecord : originalFastqReader) {
            final FastqHeader header = FastqHeader.parseFromFastqRecord(originalRecord, originalNormalizer);
            if (!header.reference().equals(refHeader)) {
                LOGGER.warn("  Invalid header in original fastq file. Record = " + originalRecord);
            } else {
                final FastqRecord recreatedMatch = recreatedFastq.get(header.key());
                if (recreatedMatch != null) {
                    if (!recreatedMatch.getReadString().equals(originalRecord.getReadString())
                            || !recreatedMatch.getBaseQualityString().equals(originalRecord.getBaseQualityString())) {
                        LOGGER.warn("  Mismatch between original and recreated fastq on record: " + originalRecord);
                        success = false;
                    }
                    recreatedFastq.remove(header.key());
                }
            }
            recordCount++;
            if (recordCount % 1E7 == 0) {
                int recordsFound = recreatedSize - recreatedFastq.size();
                LOGGER.info("  Finished mapping " + recordCount + " records. Found " + recordsFound
                        + " recreated records");
            }
        }
        final int recordsFound = recreatedSize - recreatedFastq.size();
        LOGGER.info("  Finished mapping " + recordCount + " records. Found " + recordsFound + " recreated records");

        LOGGER.info("Finished mapping records. " + recreatedFastq.size()
                + " unmapped records remaining in recreated fastq");

        return success && recreatedFastq.size() == 0;
    }

    @Nullable
    private static String referenceHeader(@NotNull final File recreatedFastqFile) {
        final FastqReader fastqReader = new FastqReader(recreatedFastqFile);
        if (fastqReader.hasNext()) {
            final FastqHeader header = FastqHeader.parseFromFastqRecord(fastqReader.next(),
                    new RecreatedFastqHeaderNormalizer());
            fastqReader.close();
            return header.reference();
        }

        return null;
    }

    @NotNull
    private static Map<FastqHeaderKey, FastqRecord> mapRecreatedFastqFile(@NotNull final File file,
            @NotNull final String refHeader, final int numRecords) {
        final FastqHeaderNormalizer normalizer = new RecreatedFastqHeaderNormalizer();
        final Map<FastqHeaderKey, FastqRecord> records = new HashMap<>(numRecords);

        final FastqReader fastqReader = new FastqReader(file);

        while (fastqReader.hasNext() && records.size() < numRecords) {
            final FastqRecord record = fastqReader.next();
            final FastqHeader header = FastqHeader.parseFromFastqRecord(record, normalizer);

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
