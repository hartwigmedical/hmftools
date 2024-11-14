package com.hartwig.hmftools.redux.write;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UMI_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.UNSET;
import static com.hartwig.hmftools.redux.write.ReadOutput.DUPLICATES;
import static com.hartwig.hmftools.redux.write.ReadOutput.MISMATCHES;
import static com.hartwig.hmftools.redux.write.ReadOutput.NONE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.FragmentStatus;

import htsjdk.samtools.SAMRecord;

public class ReadDataWriter
{
    private final ReduxConfig mConfig;
    private final BufferedWriter mWriter;

    public ReadDataWriter(final ReduxConfig config)
    {
        mConfig = config;
        mWriter = initialiseReadWriter();
    }

    public boolean enabled() { return mWriter != null; }

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

            sj.add("MapQual").add("SuppData").add("Flags");
            sj.add("FirstInPair").add("ReadReversed").add("MateReversed");
            sj.add("Unmapped").add("UnmapCoords").add("MateUnmapped").add("Supplementary").add("Secondary");

            if(mConfig.WriteReadBaseLength > 0)
            {
                sj.add("BasesStart");
                sj.add("BasesEnd");
            }

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeReadData(
            final SAMRecord read, final FragmentStatus fragmentStatus, final String fragmentCoordinates, final String umiId)
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
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(read.getReadName());
            sj.add(read.getContig());
            sj.add(String.valueOf(read.getAlignmentStart()));
            sj.add(String.valueOf(read.getAlignmentEnd()));
            sj.add(read.getCigar().toString());

            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE));
            String mateCigar = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            String unmapOrigCoords = read.getReadUnmappedFlag() ? read.getStringAttribute(UNMAP_ATTRIBUTE) : null;

            sj.add(String.valueOf(abs(read.getInferredInsertSize())));
            sj.add(read.getMateReferenceName());
            sj.add(String.valueOf(read.getMateAlignmentStart()));
            sj.add(String.valueOf(read.getDuplicateReadFlag()));
            sj.add(fragmentStatus.toString());
            sj.add(mateCigar != null ? mateCigar : "");
            sj.add(fragmentCoordinates);

            if(mConfig.UMIs.Enabled)
            {
                String umiType = read.getStringAttribute(UMI_TYPE_ATTRIBUTE);
                sj.add(umiId != null ? umiId : "");
                sj.add(umiType != null ? umiType : UmiReadType.NONE.toString());
            }

            sj.add(String.valueOf(read.getMappingQuality()));
            sj.add(suppData != null ? suppData.asDelimStr() : "N/A");
            sj.add(String.valueOf(read.getFlags()));

            boolean isPaired = read.getReadPairedFlag();
            sj.add(String.valueOf(!isPaired || read.getFirstOfPairFlag()));
            sj.add(String.valueOf(read.getReadNegativeStrandFlag()));
            sj.add(String.valueOf(isPaired && read.getMateNegativeStrandFlag()));
            sj.add(String.valueOf(read.getReadUnmappedFlag()));
            sj.add(unmapOrigCoords != null ? unmapOrigCoords : "");
            sj.add(String.valueOf(isPaired && read.getMateUnmappedFlag()));
            sj.add(String.valueOf(read.getSupplementaryAlignmentFlag()));
            sj.add(String.valueOf(read.isSecondaryAlignment()));


            if(mConfig.WriteReadBaseLength > 0 && mConfig.WriteReadBaseLength * 2 <= read.getReadBases().length)
            {
                String readBases = read.getReadString();
                sj.add(readBases.substring(0, mConfig.WriteReadBaseLength));
                sj.add(readBases.substring(readBases.length() - mConfig.WriteReadBaseLength));
            }

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            RD_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }
}
