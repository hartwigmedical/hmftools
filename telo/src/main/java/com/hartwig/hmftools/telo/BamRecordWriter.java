package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.FileWriterUtils;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamRecordWriter
{
    private final BufferedWriter mCsvWriter;
    private final SAMFileWriter mBamFileWriter;

    public BamRecordWriter(final TeloConfig config)
    {
        String sampleBamFileBasename = new File(config.SampleBamFile).getName();
        final String csvOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_read_data.csv";
        mCsvWriter = createReadDataWriter(csvOutputFile);

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.RefGenomeFile)).open(new File(config.SampleBamFile));

        final String bamOutputFile = config.OutputDir + "/" + sampleBamFileBasename + ".telo_bam.bam";
        mBamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(bamOutputFile));
    }

    public void close()
    {
        closeBufferedWriter(mCsvWriter);
        mBamFileWriter.close();
    }

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
        if(mCsvWriter != null)
        {
            try
            {
                boolean completeGroup = readGroup.isComplete();

                for(SAMRecord samRecord : readGroup.Reads)
                {
                    ReadRecord readRecord = ReadRecord.from(samRecord);

                    mCsvWriter.write(String.format("%s,%s,%d,%d,%s,%d,%s",
                            readRecord.Id, readRecord.Chromosome, readRecord.PosStart, readRecord.PosEnd,
                            readRecord.mateChromosome(), readRecord.mateStartPosition(), readRecord.hasTeloContent()));

                    String suppAlignmentData = readRecord.hasSuppAlignment() ? readRecord.getSuppAlignment().replace(",", ";") : "";
                    suppAlignmentData = suppAlignmentData.replace(",", ";");

                    mCsvWriter.write(String.format(",%s,%d,%s,%s,%s,%d",
                            readRecord.Cigar.toString(), readRecord.fragmentInsertSize(), readRecord.isFirstOfPair(),
                            readRecord.isUnmapped(), readRecord.isMateUnmapped(), readRecord.flags()));

                    mCsvWriter.write(String.format(",%s,%s,%s",
                            suppAlignmentData, completeGroup, readRecord.ReadBases));

                    mCsvWriter.newLine();
                }
            }
            catch(IOException e)
            {
                TE_LOGGER.error("failed to write read data file: {}", e.toString());
            }
        }

        if(mBamFileWriter != null)
        {
            for(SAMRecord samRecord : readGroup.Reads)
            {
                mBamFileWriter.addAlignment(samRecord);
            }
        }
    }

}
