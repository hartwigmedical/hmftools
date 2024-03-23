package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.bam.BamToolName.fromPath;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.ESVEE_FILE_ID;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.BamToolName;

import org.jetbrains.annotations.Nullable;

public final class CommonUtils
{
    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }

    public static String formOutputFile(
            final String outputDir, final String sampleId, final String fileTypeExtension, @Nullable final String outputId)
    {
        String filename = outputDir;

        filename += sampleId + "." + ESVEE_FILE_ID;

        if(outputId != null)
            filename += "." + outputId;

        filename += "." + fileTypeExtension;

        return filename;
    }

    public static void deleteInterimFile(final String filename)
    {
        try
        {
            Files.deleteIfExists(Paths.get(filename));
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error deleting interim file: {}", e.toString());
        }
    }

    public static void writeSortedBam(final String unsortedBam, final String sortedBam, final String bamToolPath, final int threads)
    {
        if(bamToolPath == null)
            return;

        SV_LOGGER.info("writing sorted BAM: {}", sortedBam);

        BamToolName toolName = fromPath(bamToolPath);

        boolean success = BamOperations.sortBam(toolName, bamToolPath, unsortedBam, sortedBam, threads);

        if(success && toolName == BamToolName.SAMTOOLS)
        {
            success = BamOperations.indexBam(toolName, bamToolPath, sortedBam, threads);

            if(success)
                deleteInterimFile(unsortedBam);
        }
    }
}
