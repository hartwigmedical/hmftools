package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.FileWriterUtils;

public class BamRecordWriter
{
    private final BufferedWriter mReadWriter;


    public BamRecordWriter(final TeloConfig config)
    {
        final String outputFile = config.OutputDir + "telo_read_data.csv";
        mReadWriter = outputFile != null ? createReadDataWriter(outputFile) : null;
    }

    public void close() { closeBufferedWriter(mReadWriter); }

    private static BufferedWriter createReadDataWriter(final String outputFile)
    {
        try
        {

            BufferedWriter writer = createBufferedWriter(outputFile, false);
            writer.write("ReadId,Chromosome,PosStart,PosEnd,MateChr,MatePosStart,HasTeloContent");
            writer.write(",Cigar,InsertSize,FirstInPair,Unmapped,MateUnmapped,Flags,SuppData,CompleteFrag,ReadBases");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            TE_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeReadGroup(final ReadGroup readGroup)
    {
        if(mReadWriter == null)
            return;

        try
        {
            boolean completeGroup = readGroup.isComplete();

            for(ReadRecord readRecord : readGroup.Reads)
            {
                mReadWriter.write(String.format("%s,%s,%d,%d,%s,%d,%s",
                        readRecord.Id, readRecord.Chromosome, readRecord.PosStart, readRecord.PosEnd,
                        readRecord.mateChromosome(), readRecord.mateStartPosition(), readRecord.hasTeloContent()));

                String suppAlignmentData = readRecord.hasSuppAlignment() ? readRecord.getSuppAlignment().replace(",", ";") : "";
                suppAlignmentData = suppAlignmentData.replace(",", ";");

                mReadWriter.write(String.format(",%s,%d,%s,%s,%s,%d",
                        readRecord.Cigar.toString(), readRecord.fragmentInsertSize(), readRecord.isFirstOfPair(),
                        readRecord.isUnmapped(), readRecord.isMateUnmapped(), readRecord.flags()));

                mReadWriter.write(String.format(",%s,%s,%s",
                        suppAlignmentData, completeGroup, readRecord.ReadBases));

                mReadWriter.newLine();
            }
        }
        catch(IOException e)
        {
            TE_LOGGER.error("failed to write read data file: {}", e.toString());
        }
    }

}
