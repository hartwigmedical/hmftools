package com.hartwig.hmftools.sullivan;

import com.google.common.collect.PeekingIterator;
import com.google.common.collect.Sets;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

public final class SullivanAlgo {

    public static void runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath) {
        Set<SullivanRecord> originalFastq = createFromFastqFile(originalFastqPath, true);
        System.out.println("Finished reading original fastq file from " + originalFastqPath +
                ". Created " + originalFastq.size() + " records.");

        Set<SullivanRecord> recreatedFastq = createFromFastqFile(recreatedFastqPath, false);
        System.out.println("Finished reading recreated fastq file from " + recreatedFastqPath +
                ". Created " + recreatedFastq.size() + " records.");

        Sets.SetView<SullivanRecord> uniqueInOriginalSet = Sets.difference(originalFastq, recreatedFastq);
        Sets.SetView<SullivanRecord> uniqueInRecreatedSet = Sets.difference(recreatedFastq, originalFastq);

        System.out.println("Found " + uniqueInOriginalSet.size() +
                " records in original that are not present in recreated");
        System.out.println("Found " + uniqueInRecreatedSet.size() +
                " records in recreated that are not present in original");
    }

    @NotNull
    private static Set<SullivanRecord> createFromFastqFile(@NotNull String path, boolean isOriginalFastq) {
        FastqHeaderParser parser = isOriginalFastq ? new OriginalFastqHeaderParser() : new RecreatedFastqHeaderParser();
        Set<SullivanRecord> records = new HashSet<SullivanRecord>();

        FastqReader fastqReader = new FastqReader(new File(path));

        for (FastqRecord record : fastqReader) {
            int recordCount = records.size();
            SullivanRecord convertedRecord = SullivanRecord.createFromFastqRecord(record, parser);
            records.add(convertedRecord);
            if (recordCount == records.size()) {
                System.out.println("WARN: Duplicate record found: " + convertedRecord);
            }
        }

        fastqReader.close();
        return records;
    }
}
