package com.hartwig.hmftools.viridian.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
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

    public static List<Integer> parseMdTagMismatchPositions(int alignmentStart, String mdTag)
    {
        List<Integer> positions = new ArrayList<>();
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
                    positions.add(position);
                }
                ++position;
            }
            else
            {
                throw new IllegalArgumentException("Invalid MD tag: " + mdTag);
            }
        }
        return positions;
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
