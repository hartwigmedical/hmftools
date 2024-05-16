package com.hartwig.hmftools.esvee.common;

import static com.hartwig.hmftools.common.bam.BamToolName.fromPath;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.PREP_FRAG_LENGTH_FILE_ID;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.BamToolName;

import org.jetbrains.annotations.Nullable;

public final class FileCommon
{
    public static final String APP_NAME = "Esvee";

    public static final String INPUT_VCF = "input_vcf";
    public static final String INPUT_VCF_DESC = "Input VCF";

    // file IDs and names
    public static final String ESVEE_FILE_ID = "esvee";
    public static final String PREP_FILE_ID = "esvee.prep";

    public static final String RAW_VCF_SUFFIX = "raw" + VCF_ZIP_EXTENSION;
    public static final String DEPTH_VCF_SUFFIX = "ref_depth" + VCF_ZIP_EXTENSION;

    public static final String FILE_NAME_DELIM = ".";

    public static final String REF_GENOME_IMAGE_EXTENSION = ".img";

    public static String formPrepInputFilename(
            final String outputDir, final String sampleId, final String fileId, @Nullable final String outputId)
    {
        return formOutputFile(outputDir, sampleId, PREP_FILE_ID, fileId, outputId);
    }

    public static String formEsveeInputFilename(
            final String outputDir, final String sampleId, final String fileId, @Nullable final String outputId)
    {
        return formOutputFile(outputDir, sampleId, ESVEE_FILE_ID, fileId, outputId);
    }

    public static String formFragmentLengthDistFilename(final String outputDir, final String sampleId)
    {
        return formPrepInputFilename(outputDir, sampleId, PREP_FRAG_LENGTH_FILE_ID, null);
    }

    public static String formOutputFile(
            final String outputDir, final String sampleId, final String appStage, final String fileType, @Nullable final String outputId)
    {
        String filename = outputDir;

        filename += sampleId + FILE_NAME_DELIM + appStage;

        if(outputId != null)
            filename += FILE_NAME_DELIM + outputId;

        filename += FILE_NAME_DELIM + fileType;

        return filename;
    }


    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
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
        }

        if(success)
            deleteInterimFile(unsortedBam);
    }
}
