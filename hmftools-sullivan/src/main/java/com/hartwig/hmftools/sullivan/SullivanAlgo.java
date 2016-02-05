package com.hartwig.hmftools.sullivan;

import com.google.common.collect.Sets;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public final class SullivanAlgo {

    private static final DateFormat DTF = new SimpleDateFormat("hh:mm:ss");

    public static boolean runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath) {
        log("Start reading original fastq file from " + originalFastqPath);
        Map<String, FastqRecord> originalFastq = createFromFastqFile(originalFastqPath, true);
        log("Finished reading original fastq file. Created " + originalFastq.size() + " records.");

        log("Start reading recreated fastq file from " + recreatedFastqPath);
        Map<String, FastqRecord> recreatedFastq = createFromFastqFile(recreatedFastqPath, false);
        log("Finished reading recreated fastq file. Created " + recreatedFastq.size() + " records.");

        boolean success = true;
        int recordCount = 0;

        for (Map.Entry<String, FastqRecord> originalMapEntry : originalFastq.entrySet()) {
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

        for (Map.Entry<String, FastqRecord> recreatedMapEntry : recreatedFastq.entrySet()) {
            log("Could not find key from recreated fastQ in original fastQ: " + recreatedMapEntry.getKey());
            success = false;
        }


        return success;
    }

    @NotNull
    private static Map<String, FastqRecord> createFromFastqFile(@NotNull String path, boolean isOriginalFastq) {
        FastqHeaderParser parser = isOriginalFastq ? new OriginalFastqHeaderParser() : new RecreatedFastqHeaderParser();
        Map<String, FastqRecord> records = new HashMap<String, FastqRecord>();

        FastqReader fastqReader = new FastqReader(new File(path));

        for (FastqRecord record : fastqReader) {
            int recordCount = records.size();
            String convertedHeader = parser.apply(record.getReadHeader());

            records.put(convertedHeader, record);
            if (recordCount == records.size()) {
                System.out.println("WARN: Duplicate record found: " + record);
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
