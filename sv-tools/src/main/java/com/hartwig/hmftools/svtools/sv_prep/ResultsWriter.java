package com.hartwig.hmftools.svtools.sv_prep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svtools.sv_prep.WriteType.BUCKET_STATS;
import static com.hartwig.hmftools.svtools.sv_prep.WriteType.READS;
import static com.hartwig.hmftools.svtools.sv_prep.WriteType.SV_BED;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

public class ResultsWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mReadWriter;
    private final BufferedWriter mBedWriter;
    private final BufferedWriter mBucketWriter;

    public ResultsWriter(final SvConfig config)
    {
        mConfig = config;

        if(mConfig.OutputDir == null)
        {
            mBucketWriter = null;
            mBedWriter = null;
            mReadWriter = null;
            return;
        }

        mBucketWriter = initialiseBucketWriter();
        mBedWriter = initialiseBedWriter();
        mReadWriter = initialiseReadWriter();
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mBedWriter);
        closeBufferedWriter(mBucketWriter);
    }

    private String formFilename(final WriteType writeType)
    {
        String filename = mConfig.OutputDir + mConfig.SampleId;

        filename += ".sv_prep.";

        if(mConfig.OutputId != null)
            filename += mConfig.OutputId + ".";

        switch(writeType)
        {
            case SV_BED: return filename + "bed";
            case READS: return filename + "reads.csv";
            case BUCKET_STATS: return filename + "buckets.csv";
            case BAM: return filename + "bam";
        }

        return null;
    }

    private BufferedWriter initialiseReadWriter()
    {
        if(!mConfig.WriteTypes.contains(READS))
            return null;

        try
        {
            String filename = formFilename(READS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("PartitionIndex,BucketId,ReadId,ReadType,GroupComplete,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",InsertSize,FragLength,MateChr,MatePosStart,FirstInPair,ReadReversed,SuppData");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeReadData(
            final ReadRecord read, int partitionIndex, int bucketId, final String readType, final boolean groupComplete)
    {
        if(mReadWriter == null)
            return;

        try
        {
            mReadWriter.write(format("%d,%d,%s,%s,%s,%s,%d,%d,%s",
                    partitionIndex, bucketId, read.Id, readType, groupComplete,
                    read.Chromosome, read.start(), read.end(), read.Cigar.toString()));

            SupplementaryReadData suppData = read.supplementaryAlignment();

            mReadWriter.write(format(",%d,%d,%s,%d,%s,%s,%s",
                    read.fragmentInsertSize(), 0, read.MateChromosome, read.MatePosStart, read.isFirstOfPair(), read.isReadReversed(),
                    suppData != null ? suppData.asCsv() : "N/A"));

            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    private BufferedWriter initialiseBedWriter()
    {
        if(!mConfig.WriteTypes.contains(SV_BED))
            return null;

        return null;
    }

    private BufferedWriter initialiseBucketWriter()
    {
        if(!mConfig.WriteTypes.contains(BUCKET_STATS))
            return null;

        try
        {
            String filename = formFilename(BUCKET_STATS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("PartitionIndex,BucketId,ReadGroups,JunctionCount,SupportingReads,InitSupportingReads,Junctions");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeBucketData(final SvBucket bucket, int partitionIndex)
    {
        if(mBucketWriter == null)
            return;

        try
        {
            StringJoiner juncDataSj = new StringJoiner(ITEM_DELIM);

            for(JunctionData junctionData : bucket.junctionPositions())
            {
                juncDataSj.add(format("%d:%d:%d:%d",
                        junctionData.Position, junctionData.Orientation, junctionData.ExactSupport, junctionData.CandidateSupport));
            }

            mBucketWriter.write(format("%d,%d,%d,%d,%d,%d,%s",
                    partitionIndex, bucket.id(), bucket.readGroups().size(), bucket.junctionPositions().size(),
                    bucket.supportingReads().size(), bucket.initialSupportingReadCount(), juncDataSj.toString()));

            mBucketWriter.newLine();;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }


}
