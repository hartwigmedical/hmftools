package com.hartwig.hmftools.markdups.write;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.markdups.common.FragmentStatus.UNSET;
import static com.hartwig.hmftools.markdups.write.ReadOutput.DUPLICATES;
import static com.hartwig.hmftools.markdups.write.ReadOutput.MISMATCHES;
import static com.hartwig.hmftools.markdups.write.ReadOutput.NONE;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.markdups.MarkDupsConfig;
import com.hartwig.hmftools.markdups.common.FragmentStatus;

import htsjdk.samtools.SAMRecord;

public class ReadDataWriter
{
    private final MarkDupsConfig mConfig;
    private final BufferedWriter mWriter;

    public ReadDataWriter(final MarkDupsConfig config)
    {
        mConfig = config;
        mWriter = initialiseReadWriter();
    }

    private BufferedWriter initialiseReadWriter()
    {
        if(mConfig.LogReadType == NONE)
            return null;

        try
        {
            String filename = mConfig.formFilename("reads");
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId\tChromosome\tPosStart\tPosEnd\tCigar");
            writer.write("\tInsertSize\tMateChr\tMatePosStart\tDuplicate\tCalcDuplicate\tMateCigar\tCoords");

            if(mConfig.UMIs.Enabled)
                writer.write("\tUmi\tUmiType");

            writer.write("\tAvgBaseQual\tMapQual\tSuppData\tFlags");
            writer.write("\tFirstInPair\tReadReversed\tUnmapped\tMateUnmapped\tSupplementary\tSecondary");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeReadData(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates,
            final double avgBaseQual, final String umiId)
    {
        if(mWriter == null)
            return;

        if(mConfig.LogReadType == DUPLICATES)
        {
            if(!read.getDuplicateReadFlag() && !fragmentStatus.isDuplicate())
                return;
        }
        else if(mConfig.LogReadType == MISMATCHES)
        {
            if(fragmentStatus != UNSET)
            {
                if(read.getDuplicateReadFlag() == (fragmentStatus == DUPLICATE))
                    return;
            }
        }

        try
        {
            mWriter.write(format("%s\t%s\t%d\t%d\t%s",
                    read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(), read.getCigar()));

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));

            mWriter.write(format("\t%d\t%s\t%d\t%s\t%s\t%s\t%s",
                    abs(read.getInferredInsertSize()), read.getMateReferenceName(), read.getMateAlignmentStart(),
                    read.getDuplicateReadFlag(), fragmentStatus, read.hasAttribute(MATE_CIGAR_ATTRIBUTE), fragmentCoordinates));

            if(mConfig.UMIs.Enabled)
            {
                String umiType = read.getStringAttribute(UMI_TYPE_ATTRIBUTE);
                mWriter.write(format("\t%s\t%s", umiId != null ? umiId : "", umiType != null ? umiType : UmiReadType.NONE));
            }

            mWriter.write(format("\t%.2f\t%d\t%s\t%d",
                    avgBaseQual, read.getMappingQuality(), suppData != null ? suppData.asCsv() : "N/A", read.getFlags()));

            mWriter.write(format("\t%s\t%s\t%s\t%s\t%s\t%s",
                    read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), read.getReadUnmappedFlag(),
                    read.getMateUnmappedFlag(), read.getSupplementaryAlignmentFlag(), read.isSecondaryAlignment()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            MD_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }
}
