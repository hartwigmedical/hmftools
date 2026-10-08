package com.hartwig.hmftools.viridian.common;

import static java.nio.file.StandardCopyOption.ATOMIC_MOVE;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

public final class AtomicOutput
{
    @FunctionalInterface
    public interface Writer<T>
    {
        T writeTo(String file) throws IOException;
    }

    private static final String TEMP_FILE_PREFIX = "tmp.";

    // Deletes any previous file, writes to a temporary sibling, and moves it into place if successful.
    public static <T> T write(String file, Writer<T> writer) throws IOException
    {
        Path finalFile = Path.of(file);
        Path tempFile = tempFile(finalFile);

        Files.deleteIfExists(finalFile);
        T result = writeOrDiscard(writer, tempFile, List.of(tempFile));
        Files.move(tempFile, finalFile, ATOMIC_MOVE);
        return result;
    }

    // Writes a BAM and its index atomically. The writer must write the BAM with htsjdk index creation enabled, which names
    // the index after the BAM's stem; it is moved to the conventional "<bam>.bai". The index moves first, so the BAM
    // existing means both are complete.
    // The previous BAM and index are deleted too.
    public static <T> T writeIndexedBam(String bamFile, Writer<T> writer) throws IOException
    {
        if(!bamFile.endsWith(BAM_EXTENSION))
        {
            throw new IllegalArgumentException("Not a BAM file name: " + bamFile);
        }

        Path finalBam = Path.of(bamFile);
        Path finalIndex = Path.of(bamFile + BAM_INDEX_EXTENSION);
        Path tempBam = tempFile(finalBam);
        Path tempIndex = htsjdkIndexFile(tempBam);

        Files.deleteIfExists(finalIndex);
        Files.deleteIfExists(finalBam);
        T result = writeOrDiscard(writer, tempBam, List.of(tempBam, tempIndex));
        Files.move(tempIndex, finalIndex, ATOMIC_MOVE);
        Files.move(tempBam, finalBam, ATOMIC_MOVE);
        return result;
    }

    private static Path tempFile(Path file)
    {
        return file.resolveSibling(TEMP_FILE_PREFIX + file.getFileName());
    }

    private static Path htsjdkIndexFile(Path bamFile)
    {
        String fileName = bamFile.getFileName().toString();
        String stem = fileName.substring(0, fileName.length() - BAM_EXTENSION.length());
        return bamFile.resolveSibling(stem + BAM_INDEX_EXTENSION);
    }

    private static <T> T writeOrDiscard(Writer<T> writer, Path tempFile, List<Path> tempFiles) throws IOException
    {
        try
        {
            return writer.writeTo(tempFile.toString());
        }
        catch(Throwable writeFailure)
        {
            for(Path file : tempFiles)
            {
                try
                {
                    Files.deleteIfExists(file);
                }
                catch(IOException deleteFailure)
                {
                    writeFailure.addSuppressed(deleteFailure);
                }
            }
            throw writeFailure;
        }
    }
}
