package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.BamToolName.fromPath;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateUnmapped;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.common.SvConstants.FILE_NAME_DELIM;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.BamToolName;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterConfig;
import com.hartwig.hmftools.esvee.read.Read;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public final class CommonUtils
{
    public static boolean isDiscordantFragment(
            final SAMRecord read, final int fragmentLengthUpperBound, @Nullable final SupplementaryReadData suppData)
    {
        if(read.getReadUnmappedFlag() || !read.getReadPairedFlag() || read.getMateUnmappedFlag())
            return false;

        // supplementaries need to check their primary read chromosomes, not their own
        if(read.getSupplementaryAlignmentFlag() && suppData != null)
        {
            if(!suppData.Chromosome.equals(read.getMateReferenceName()))
                return true;
        }
        else if(!read.getReferenceName().equals(read.getMateReferenceName()))
        {
            return true;
        }

        // inversion
        if(read.getReadNegativeStrandFlag() == read.getMateNegativeStrandFlag())
            return true;

        int fragmentSize = abs(read.getInferredInsertSize());

        return fragmentSize == 0 || fragmentSize >= fragmentLengthUpperBound;
    }

    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public static BamSlicer createBamSlicer()
    {
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
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

    public static void copyArray(final byte[] source, final byte[] dest, final int sourceIndexStart, final int destIndexStart)
    {
        int d = destIndexStart;
        for(int s = sourceIndexStart; s < source.length && d < dest.length; ++s, ++d)
        {
            dest[d] = source[s];
        }
    }

    public static byte[] copyArray(final byte[] source)
    {
        byte[] dest = new byte[source.length];

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }

        return dest;
    }

    public static byte[] subsetArray(final byte[] source, final int startIndex, final int endIndex)
    {
        byte[] dest = new byte[endIndex - startIndex + 1];

        int newIndex = 0;
        for(int index = startIndex; index <= endIndex; ++index, ++newIndex)
        {
            dest[newIndex] = source[index];
        }

        return dest;
    }

    public static int[] copyArray(final int[] source)
    {
        int[] dest = new int[source.length];

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }

        return dest;
    }

    public static byte[] addByteArray(final byte[] first, final byte[] second)
    {
        byte[] combined = new byte[first.length + second.length];

        for(int i = 0; i < first.length; ++i)
        {
            combined[i] = first[i];
        }

        for(int i = 0; i < second.length; ++i)
        {
            combined[first.length + i] = second[i];
        }

        return combined;
    }

    public static byte[] reverseBytes(final byte[] bases)
    {
        String reversed = Nucleotides.reverseComplementBases(new String(bases));
        return reversed.getBytes();
    }

    public static void initialise(final byte[] array, final byte value)
    {
        for(int i = 0; i < array.length; ++i)
        {
            array[i] = value;
        }
    }

    public static void initialise(final int[] array, final int value)
    {
        for(int i = 0; i < array.length; ++i)
        {
            array[i] = value;
        }
    }

}
