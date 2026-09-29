package com.hartwig.hmftools.viridian.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

public class Utils
{
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
