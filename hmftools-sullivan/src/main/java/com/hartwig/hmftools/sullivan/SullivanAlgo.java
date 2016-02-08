package com.hartwig.hmftools.sullivan;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public final class SullivanAlgo {

    private static final DateFormat DTF = new SimpleDateFormat("hh:mm:ss");

    public static boolean runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath) {
        return runSullivanAlgo(originalFastqPath, recreatedFastqPath, Integer.MAX_VALUE);
    }

    public static boolean runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath,
                                          int numRecords) {
        String refHeader = referenceHeader(originalFastqPath);
        if (refHeader == null) {
            log("No ref header could be isolated from fastq file on " + originalFastqPath);
            return false;
        } else {
            log("Generated ref header: " + refHeader);
        }

        log("Start reading original fastq file from " + originalFastqPath);
        Map<FastqHeaderKey, FastqRecord> originalFastq = mapOriginalFastqFile(originalFastqPath, refHeader, numRecords);
        log("Finished reading original fastq file. Created " + originalFastq.size() + " records.");

        FastqReader recreatedFastqReader = createFastqReader(recreatedFastqPath);
        FastqHeaderNormalizer recreatedNormalizer = new RecreatedFastqHeaderNormalizer();

        boolean success = true;
        int recordCount = 0;

        log("Start mapping process from recreated to original fastq");
        for (FastqRecord recreatedRecord : recreatedFastqReader) {
            FastqHeader header = FastqHeader.parseFromFastqRecord(recreatedRecord, recreatedNormalizer);
            if (!header.reference().equals(refHeader)) {
                log("  Invalid header in recreated fastq file. Record = " + recreatedRecord);
            } else {
                FastqRecord originalMatch = originalFastq.get(header.key());
                if (originalMatch != null) {
                    if (!originalMatch.getReadString().equals(recreatedRecord.getReadString()) ||
                        !originalMatch.getBaseQualityString().equals(recreatedRecord.getBaseQualityString())) {
                        log("  Mismatch between original and recreated fastq on record: " + recreatedRecord);
                        success = false;
                    }
                    originalFastq.remove(header.key());
                }
            }
            recordCount++;
            if (recordCount % 1E6 == 0) {
                log("  Finished mapping " + recordCount + " records.");
            }
        }

        log("Finished mapping records. " + originalFastq.size() + " unmapped records remaining in original fastq");

        return success && originalFastq.size() == 0;
    }

    @Nullable
    private static String referenceHeader(@NotNull String originalFastqPath) {
        FastqReader fastqReader = createFastqReader(originalFastqPath);
        if (fastqReader.hasNext()) {
            FastqHeader header = FastqHeader.parseFromFastqRecord(fastqReader.next(), new OriginalFastqHeaderNormalizer());
            fastqReader.close();
            return header.reference();
        }

        return null;
    }

    @NotNull
    private static Map<FastqHeaderKey, FastqRecord> mapOriginalFastqFile(
            @NotNull String path, @NotNull String refHeader, int numRecords) {
        FastqHeaderNormalizer normalizer = new OriginalFastqHeaderNormalizer();
        Map<FastqHeaderKey, FastqRecord> records = new HashMap<FastqHeaderKey, FastqRecord>();

        FastqReader fastqReader = createFastqReader(path);

        while (fastqReader.hasNext() && records.size() < numRecords) {
            FastqRecord record = fastqReader.next();
            FastqHeader header = FastqHeader.parseFromFastqRecord(record, normalizer);

            if (records.containsKey(header.key())) {
                log("  WARN: Duplicate record found: " + record);
            }

            records.put(header.key(), record);

            if (!header.reference().equals(refHeader)) {
                log("  WARN: Header mismatch with reference header: " + header.reference());
            }

            if (records.size() % 1E6 == 0) {
                log("  Finished reading " + records.size() + " records");
            }
        }

        fastqReader.close();
        return records;
    }

    @NotNull
    private static FastqReader createFastqReader(@NotNull String path) {
        return new FastqReader(new File(path));
    }

    private static void log(@NotNull String info) {
        System.out.println(DTF.format(new Date()) + ": " + info);
    }
}
