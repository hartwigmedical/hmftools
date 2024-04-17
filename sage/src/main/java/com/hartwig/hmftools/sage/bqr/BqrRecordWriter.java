package com.hartwig.hmftools.sage.bqr;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.samtools.SAMRecord;

public class BqrRecordWriter
{
    private final SageConfig mConfig;
    private final BufferedWriter mReadWriter;
    private final BufferedWriter mPositionWriter;

    public BqrRecordWriter(final SageConfig config, final String sampleId)
    {
        mConfig = config;
        mPositionWriter = initialisePositionWriter(sampleId);
        mReadWriter = initialiseReadWriter(sampleId);
    }

    public boolean enabled() { return mReadWriter != null || mPositionWriter != null; }

    private BufferedWriter initialisePositionWriter(final String sampleId)
    {
        if(!mConfig.BQR.WritePositions)
            return null;

        try
        {
            String outputFile = mConfig.outputDir() + sampleId + ".bqr_positions.tsv";;
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Chromosome").add("Position").add("Ref").add("Alt").add("TriNuc");
            sj.add("Quality").add("ReadType").add("Count");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            e.printStackTrace();
            SG_LOGGER.error("failed to write BQR position data", e.toString());
            System.exit(1);
        }

        return null;
    }
    public synchronized void writePositionData(
            final String chromosome, final int position, final byte ref, final byte alt,
            final byte[] trinucleotideContext, final byte quality, final BqrReadType readType, final int count)
    {
        if(mPositionWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(chromosome);
            sj.add(String.valueOf(position));
            sj.add(String.valueOf((char)ref));
            sj.add(String.valueOf((char)alt));
            sj.add(new String(trinucleotideContext));
            sj.add(String.valueOf(quality));
            sj.add(String.valueOf(readType));
            sj.add(String.valueOf(count));

            mPositionWriter.write(sj.toString());
            mPositionWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write BQR position data", e.toString());
        }
    }

    private BufferedWriter initialiseReadWriter(final String sampleId)
    {
        if(!mConfig.BQR.WriteReads)
            return null;

        try
        {
            String outputFile = mConfig.outputDir() + sampleId + ".bqr_reads.tsv";;
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Ref").add("Alt").add("TriNuc");
            sj.add("Chromosome").add("Position").add("ReadIndex").add("Quality").add("ReadType");
            sj.add("ReadOrient").add("FragOrient").add("FivePrimeDist").add("ThreePrimeDist");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            e.printStackTrace();
            SG_LOGGER.error("failed to write BQR position data", e.toString());
            System.exit(1);
        }

        return null;
    }

    public synchronized void writeRecordData(
            final SAMRecord record, final int position, final int readIndex, final byte ref, final byte alt,
            final byte[] trinucleotideContext, final byte quality, final BqrReadType readType)
    {
        if(mReadWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(String.valueOf((char)ref));
            sj.add(String.valueOf((char)alt));
            sj.add(new String(trinucleotideContext));
            sj.add(record.getReferenceName());
            sj.add(String.valueOf(position));
            sj.add(String.valueOf(readIndex));
            sj.add(String.valueOf(quality));
            sj.add(String.valueOf(readType));

            boolean isNegStrand = record.getReadNegativeStrandFlag();
            int readLength = record.getReadBases().length;
            int distanceFromFivePrimeEnd = isNegStrand ? readLength - readIndex - 1 : readIndex;
            int distanceFromThreePrimeEnd = readLength - distanceFromFivePrimeEnd;

            sj.add(String.valueOf(orientationChar(record.getReadNegativeStrandFlag())));
            sj.add(fragmentOrientationStr(record));
            sj.add(String.valueOf(distanceFromFivePrimeEnd));
            sj.add(String.valueOf(distanceFromThreePrimeEnd));

            mReadWriter.write(sj.toString());
            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write BQR read data", e.toString());
        }
    }

    private static String fragmentOrientationStr(final SAMRecord record)
    {
        if(!record.getReadPairedFlag())
            return orientationChar(record.getReadNegativeStrandFlag()) + "1";

        if(record.getFirstOfPairFlag())
            return format("%c1%c2", orientationChar(record.getReadNegativeStrandFlag()), orientationChar(record.getMateNegativeStrandFlag()));
        else
            return format("%c1%c2", orientationChar(record.getMateNegativeStrandFlag()), orientationChar(record.getReadNegativeStrandFlag()));
    }

    private static char orientationChar(boolean negStrand) { return negStrand ? 'R' : 'F'; }

    public void close()
    {
        closeBufferedWriter(mPositionWriter);
        closeBufferedWriter(mReadWriter);
    }
}
