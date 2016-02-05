package com.hartwig.hmftools.sullivan;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public final class SullivanAlgo {

    private static final DateFormat DTF = new SimpleDateFormat("hh:mm:ss");
    private static final FastqHeaderNormalizer ORIGINAL_HEADER_NORMALIZER = new OriginalFastqHeaderNormalizer();
    private static final FastqHeaderNormalizer RECREATED_HEADER_NORMALIZER = new RecreatedFastqHeaderNormalizer();

    public static boolean runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath) {
        String refHeader = referenceHeader(originalFastqPath);
        log("Generated ref header: " + refHeader);

        log("Start reading original fastq file from " + originalFastqPath);
        Map<FastqHeaderKey, FastqRecord> originalFastq = createFromFastqFile(originalFastqPath, refHeader, true);
        log("Finished reading original fastq file. Created " + originalFastq.size() + " records.");

        log("Start reading recreated fastq file from " + recreatedFastqPath);
        Map<FastqHeaderKey, FastqRecord> recreatedFastq = createFromFastqFile(recreatedFastqPath, refHeader, false);
        log("Finished reading recreated fastq file. Created " + recreatedFastq.size() + " records.");

        boolean success = true;
        int recordCount = 0;

        for (Map.Entry<FastqHeaderKey, FastqRecord> originalMapEntry : originalFastq.entrySet()) {
            FastqRecord originalEntry = originalMapEntry.getValue();
            FastqRecord recreatedEntry = recreatedFastq.get(originalMapEntry.getKey());
            if (recreatedEntry == null) {
                log("  Could not find following original entry in recreated fastq: " + originalMapEntry.getKey());
                success = false;
            } else if (!recreatedEntry.getReadString().equals(originalEntry.getReadString()) ||
                    !recreatedEntry.getBaseQualityString().equals(originalEntry.getBaseQualityString())) {
                log("  Mismatch in entry for header " + originalMapEntry.getKey());
                success = false;
            }

            if (recreatedEntry != null) {
                recreatedFastq.remove(originalMapEntry.getKey());
            }

            recordCount++;
            if (recordCount % 1E6 == 0) {
                log("  Finished comparing " + recordCount + " records.");
            }
        }

        for (Map.Entry<FastqHeaderKey, FastqRecord> recreatedMapEntry : recreatedFastq.entrySet()) {
            log("Could not find key from recreated fastQ in original fastQ: " +
                    recreatedMapEntry.getValue().getReadHeader());
            success = false;
        }

        return success;
    }

    @NotNull
    private static String referenceHeader(@NotNull String originalFastqPath) {
        FastqReader fastqReader = new FastqReader(new File(originalFastqPath));
        FastqHeader header = FastqHeader.parseFromFastqRecord(fastqReader.next(), new OriginalFastqHeaderNormalizer());
        fastqReader.close();
        return header.reference();
    }

    @NotNull
    private static Map<FastqHeaderKey, FastqRecord> createFromFastqFile(
            @NotNull String path, @NotNull String refHeader, boolean isOriginalFastq) {
        FastqHeaderNormalizer normalizer =
                isOriginalFastq ? new OriginalFastqHeaderNormalizer() : new RecreatedFastqHeaderNormalizer();
        Map<FastqHeaderKey, FastqRecord> records = new HashMap<FastqHeaderKey, FastqRecord>();

        FastqReader fastqReader = new FastqReader(new File(path));

        for (FastqRecord record : fastqReader) {
            int recordCount = records.size();
            FastqHeader header = FastqHeader.parseFromFastqRecord(record, normalizer);

            records.put(header.key(), record);
            if (recordCount == records.size()) {
                log("  WARN: Duplicate record found: " + record);
            }

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

    private static void log(@NotNull String info) {
        System.out.println(DTF.format(new Date()) + ": " + info);
    }
}
