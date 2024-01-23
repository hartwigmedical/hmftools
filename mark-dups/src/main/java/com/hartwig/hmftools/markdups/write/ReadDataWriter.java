package com.hartwig.hmftools.markdups.write;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
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
import java.util.StringJoiner;

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

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("ReadId").add("Chromosome").add("PosStart").add("PosEnd").add("Cigar");
            sj.add("InsertSize").add("MateChr").add("MatePosStart").add("Duplicate").add("CalcDuplicate").add("MateCigar").add("Coords");

            if(mConfig.UMIs.Enabled)
                sj.add("Umi").add("UmiType");

            sj.add("AvgBaseQual").add("MapQual").add("SuppData").add("Flags");
            sj.add("FirstInPair").add("ReadReversed").add("MateReversed");
            sj.add("Unmapped").add("UnmapCoords").add("MateUnmapped").add("Supplementary").add("Secondary");

            writer.write(sj.toString());
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
            String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            String unmapOrigCoords = read.getReadUnmappedFlag() ? read.getStringAttribute(UNMAP_ATTRIBUTE) : null;

            mWriter.write(format("\t%d\t%s\t%d\t%s\t%s\t%s\t%s",
                    abs(read.getInferredInsertSize()), read.getMateReferenceName(), read.getMateAlignmentStart(),
                    read.getDuplicateReadFlag(), fragmentStatus, mateCigar != null ? mateCigar : "", fragmentCoordinates));

            if(mConfig.UMIs.Enabled)
            {
                String umiType = read.getStringAttribute(UMI_TYPE_ATTRIBUTE);
                mWriter.write(format("\t%s\t%s", umiId != null ? umiId : "", umiType != null ? umiType : UmiReadType.NONE));
            }

            mWriter.write(format("\t%.2f\t%d\t%s\t%d",
                    avgBaseQual, read.getMappingQuality(), suppData != null ? suppData.asCsv() : "N/A", read.getFlags()));

            boolean isPaired = read.getReadPairedFlag();
            mWriter.write(format("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    !isPaired || read.getFirstOfPairFlag(), read.getReadNegativeStrandFlag(), isPaired && read.getMateNegativeStrandFlag(),
                    read.getReadUnmappedFlag(), unmapOrigCoords != null ? unmapOrigCoords : "",
                    isPaired && read.getMateUnmappedFlag(), read.getSupplementaryAlignmentFlag(), read.isSecondaryAlignment()));

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
