package com.hartwig.hmftools.sullivan;

import com.google.common.annotations.VisibleForTesting;
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

    public static boolean runSullivanAlgo(@NotNull String originalFastqPath, @NotNull String recreatedFastqPath,
                                          @Nullable String mergeOrigFastqPath, boolean isDirectoryMode, int numRecords) {
        File originalPath = new File(originalFastqPath);
        File recreatedPath = new File(recreatedFastqPath);

        File[] originalFiles;
        File[] recreatedFiles;
        if (isDirectoryMode) {
            log("Running sullivan algo in directory mode");
            assert originalPath.isDirectory() && recreatedPath.isDirectory();
            originalFiles = originalPath.listFiles();
            assert originalFiles != null;

            File[] mergeOrigFiles = new File[]{};
            if (mergeOrigFastqPath != null) {
                mergeOrigFiles = (new File(mergeOrigFastqPath)).listFiles();
                assert mergeOrigFiles != null;
            }

            recreatedFiles = new File[originalFiles.length + mergeOrigFiles.length];
            for (int i = 0; i < originalFiles.length; i++) {
                recreatedFiles[i] = new File(recreatedFastqPath + File.separator +
                        fromOriginalToRecreatedFileName(originalFiles[i].getName()));
                log("Mapped " + originalFiles[i].getPath() + " to " + recreatedFiles[i].getPath());
            }

            for (int i = 0; i < mergeOrigFiles.length; i++) {
                recreatedFiles[i + originalFiles.length] = new File(recreatedFastqPath + File.separator +
                        fromOriginalToRecreatedFileName(mergeOrigFiles[i].getName()));
                log("Mapped " + mergeOrigFiles[i].getPath() + " to " +
                        recreatedFiles[i + originalFiles.length].getPath());
            }
        } else {
            log("Running sullivan algo in file mode");
            assert originalPath.isFile() && recreatedPath.isFile();
            originalFiles = new File[]{originalPath};
            recreatedFiles = new File[]{recreatedPath};
        }

        assert originalFiles.length == recreatedFiles.length;

        boolean success = true;
        for (int i = 0; i < originalFiles.length; i++) {
            success = success &&
                    runSullivanOnTwoFiles(originalFiles[i], recreatedFiles[i], numRecords);
        }
        return success;
    }

    @VisibleForTesting
    @NotNull
    static String fromOriginalToRecreatedFileName(@NotNull String name) {
        String nameWithoutExtension = name.substring(0, name.indexOf("."));

        String splitRegExp = "_";
        String[] parts = nameWithoutExtension.split(splitRegExp);
        String readGroup = parts[4].equals("R1") ? "1" : "2";
        return parts[0] + splitRegExp + parts[1] + splitRegExp + parts[2] + splitRegExp + parts[3] +
                splitRegExp + parts[5] + splitRegExp + readGroup + ".fastq";
    }

    private static boolean runSullivanOnTwoFiles(@NotNull File originalFastqFile, @NotNull File recreatedFastqFile,
                                         int numRecords) {
        assert originalFastqFile.isFile() && recreatedFastqFile.isFile();

        log("Start reading original fastq file from " + originalFastqFile.getPath());

        String refHeader = referenceHeader(originalFastqFile);
        if (refHeader == null) {
            log("No ref header could be isolated from fastq file on " + originalFastqFile.getName());
            return false;
        } else {
            log("Generated ref header: " + refHeader);
        }

        Map<FastqHeaderKey, FastqRecord> originalFastq = mapOriginalFastqFile(originalFastqFile, refHeader, numRecords);
        int originalSize = originalFastq.size();
        log("Finished reading original fastq file. Created " + originalSize + " records.");

        FastqReader recreatedFastqReader = new FastqReader(recreatedFastqFile);
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
            if (recordCount % 1E7 == 0) {
                int recordsFound = originalSize - originalFastq.size();
                log("  Finished mapping " + recordCount + " records. Found " + recordsFound + " original records");
            }
        }
        int recordsFound = originalSize - originalFastq.size();
        log("  Finished mapping " + recordCount + " records. Found " + recordsFound + " original records");

        log("Finished mapping records. " + originalFastq.size() + " unmapped records remaining in original fastq");

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
                log("  WARN: Duplicate record found: " + record);
            }

            records.put(header.key(), record);

            if (!header.reference().equals(refHeader)) {
                log("  WARN: Header mismatch with reference header: " + header.reference());
            }

            if (records.size() % 1E7 == 0) {
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
