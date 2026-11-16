package com.hartwig.hmftools.viridian.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.function.BiFunction;

import htsjdk.samtools.SAMRecord;

public class Utils
{
    public static <T> T getRequiredSamTag(SAMRecord record, String tag, BiFunction<SAMRecord, String, T> getter)
    {
        T value = getter.apply(record, tag);
        if(value == null)
        {
            throw new IllegalStateException(
                    String.format("aligned read missing %s tag on contig: %s", tag, record.getReferenceName()));
        }
        return value;
    }

    // Cached value to save memory on empty arrays which are common.
    private static final int[] NO_POSITIONS = new int[0];

    public static int[] parseMdTagMismatchPositions(int alignmentStart, String mdTag)
    {
        int[] positions = new int[mdTag.length()];
        int positionCount = 0;
        int position = alignmentStart;
        int matchedBases = 0;
        boolean inDeletion = false;
        for(int i = 0; i < mdTag.length(); ++i)
        {
            char c = mdTag.charAt(i);
            if(c >= '0' && c <= '9')
            {
                matchedBases = matchedBases * 10 + (c - '0');
                inDeletion = false;
                continue;
            }

            position += matchedBases;
            matchedBases = 0;
            if(c == '^')
            {
                inDeletion = true;
            }
            else if(Character.isLetter(c))
            {
                if(!inDeletion)
                {
                    positions[positionCount++] = position;
                }
                ++position;
            }
            else
            {
                throw new IllegalArgumentException("Invalid MD tag: " + mdTag);
            }
        }
        return positionCount == 0 ? NO_POSITIONS : Arrays.copyOf(positions, positionCount);
    }

    // htsjdk names the index after the BAM's file stem, but tooling looks for it appended to the full name.
    // So "HG002.viridian.representative.bai" is renamed to "HG002.viridian.representative.bam.bai".
    public static void fixBamIndexName(String bamFile)
    {
        Path written = Path.of(bamFile.substring(0, bamFile.length() - BAM_EXTENSION.length()) + BAM_INDEX_EXTENSION);
        try
        {
            Files.move(written, Path.of(bamFile + BAM_INDEX_EXTENSION), StandardCopyOption.REPLACE_EXISTING);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to rename BAM index: " + written, e);
        }
    }
}
