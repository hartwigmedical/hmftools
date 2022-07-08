package com.hartwig.hmftools.svprep;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BUCKET_STATS;
import static com.hartwig.hmftools.svprep.WriteType.JUNCTIONS;
import static com.hartwig.hmftools.svprep.WriteType.READS;
import static com.hartwig.hmftools.svprep.WriteType.SV_BED;

import static htsjdk.samtools.SAMFlag.MATE_UNMAPPED;
import static htsjdk.samtools.SAMFlag.PROPER_PAIR;
import static htsjdk.samtools.SAMFlag.READ_UNMAPPED;
import static htsjdk.samtools.SAMFlag.SECONDARY_ALIGNMENT;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.svprep.reads.JunctionData;
import com.hartwig.hmftools.svprep.reads.ReadFilterType;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.RemoteJunction;
import com.hartwig.hmftools.svprep.reads.BucketData;

public class ResultsWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mReadWriter;
    private final BufferedWriter mBedWriter;
    private final BufferedWriter mBucketWriter;
    private final BufferedWriter mJunctionWriter;
    private final BamWriter mBamWriter;

    public ResultsWriter(final SvConfig config)
    {
        mConfig = config;

        if(mConfig.OutputDir == null)
        {
            mBucketWriter = null;
            mBedWriter = null;
            mReadWriter = null;
            mJunctionWriter = null;
            mBamWriter = null;
            return;
        }

        mBucketWriter = initialiseBucketWriter();
        mJunctionWriter = initialiseJunctionWriter();
        mBedWriter = initialiseBedWriter();
        mReadWriter = initialiseReadWriter();
        mBamWriter = new BamWriter(config);
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mBedWriter);
        closeBufferedWriter(mBucketWriter);
        closeBufferedWriter(mJunctionWriter);
        mBamWriter.close();
    }

    private BufferedWriter initialiseReadWriter()
    {
        if(!mConfig.WriteTypes.contains(READS))
            return null;

        try
        {
            String filename = mConfig.formFilename(READS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("PartitionIndex,JunctionPosition,ReadId,ReadType,GroupComplete,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",FragLength,MateChr,MatePosStart,FirstInPair,ReadReversed,SuppData");
            writer.write(",MapQual,MultiMapped,Flags,Proper,Duplicate,Unmapped,MateUnmapped,Secondary,Supplementary");

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
            final ReadRecord read, int partitionIndex, int junctionPosition, final String readType, final boolean groupComplete)
    {
        if(mReadWriter == null)
            return;

        try
        {
            mReadWriter.write(format("%d,%d,%s,%s,%s,%s,%d,%d,%s",
                    partitionIndex, junctionPosition, read.id(), readType, groupComplete,
                    read.Chromosome, read.start(), read.end(), read.cigar().toString()));

            SupplementaryReadData suppData = read.supplementaryAlignment();

            mReadWriter.write(format(",%d,%s,%d,%s,%s,%s",
                    read.fragmentInsertSize(), read.MateChromosome, read.MatePosStart, read.isFirstOfPair(), read.isReadReversed(),
                    suppData != null ? suppData.asCsv() : "N/A"));

            mReadWriter.write(format(",%d,%s,%d,%s,%s,%s,%s,%s,%s",
                    read.MapQuality, read.isMultiMapped(), read.flags(),
                    read.hasFlag(PROPER_PAIR), read.isDuplicate(), read.hasFlag(READ_UNMAPPED), read.hasFlag(MATE_UNMAPPED),
                    read.hasFlag(SECONDARY_ALIGNMENT), read.hasFlag(SUPPLEMENTARY_ALIGNMENT)));

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
            String filename = mConfig.formFilename(BUCKET_STATS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,PosStart,PosEnd,PartitionIndex,BucketId");
            writer.write(",ReadGroups,CompleteGroups,JunctionCount,SupportingReads,InitSupportingReads,Junctions");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create bucket writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeBucketData(final BucketData bucket, int partitionIndex)
    {
        if(mBucketWriter == null)
            return;

        try
        {
            StringJoiner juncDataSj = new StringJoiner(ITEM_DELIM);

            Collections.sort(bucket.junctions(), new JunctionData.JunctionDataSorter());

            for(int i = 0; i < min(bucket.junctions().size(), 20); ++i)
            {
                JunctionData junctionData = bucket.junctions().get(i);
                juncDataSj.add(format("%d:%d:%d:%d",
                        junctionData.Position, junctionData.Orientation, junctionData.exactFragmentCount(), junctionData.supportingReadCount()));
            }

            mBucketWriter.write(format("%s,%d,%d,%d,%d",
                    bucket.region().Chromosome, bucket.region().start(), bucket.region().end(), partitionIndex, bucket.id()));

            long completeGroups = bucket.readGroups().stream().filter(x -> x.isComplete()).count();
            mBucketWriter.write(format(",%d,%d,%d,%d,%d,%s",
                    bucket.readGroups().size(), completeGroups, bucket.junctions().size(),
                    bucket.supportingReads().size(), bucket.initialSupportingReadCount(), juncDataSj.toString()));

            mBucketWriter.newLine();;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write bucket data: {}", e.toString());
        }
    }

    private BufferedWriter initialiseJunctionWriter()
    {
        if(!mConfig.WriteTypes.contains(BUCKET_STATS))
            return null;

        try
        {
            String filename = mConfig.formFilename(JUNCTIONS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,BucketStart,BucketEnd,Position,Orientation,Fragments,SupportingReads,LowMapQualFrags,Hotspot,InitialReadId");
            writer.write(",RemoteJunctionCount,RemoteChromosome,RemotePosition,RemoteOrientation");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create junction writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeJunctionData(final String chromosome, final List<JunctionData> junctions)
    {
        if(mJunctionWriter == null)
            return;

        try
        {
            for(JunctionData junctionData : junctions)
            {
                int lowMapQualFrags = (int) junctionData.JunctionGroups.stream()
                        .filter(x -> x.reads().stream().anyMatch(y -> y.filters() == ReadFilterType.MIN_MAP_QUAL.flag())).count();

                String junctionStr = format("%s,%d,%d,%d,%d,%d,%s,%s,%d",
                        chromosome, junctionData.Position,
                        junctionData.Orientation, junctionData.exactFragmentCount(), junctionData.supportingReadCount(), lowMapQualFrags,
                        junctionData.hotspot(), junctionData.InitialRead.id(), junctionData.RemoteJunctions.size());

                if(!junctionData.RemoteJunctions.isEmpty())
                {
                    for(RemoteJunction remoteJunction : junctionData.RemoteJunctions)
                    {
                        mJunctionWriter.write(format("%s,%s,%d,%d",
                                junctionStr, remoteJunction.Chromosome, remoteJunction.Position, remoteJunction.Orientation));
                        mJunctionWriter.newLine();
                        ;
                    }
                }
                else
                {
                    mJunctionWriter.write(format("%s,,,", junctionStr));
                    mJunctionWriter.newLine();
                    ;
                }
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write junction data: {}", e.toString());
        }
    }

    public synchronized void writeBamRecords(final BucketData bucket)
    {
        if(mBamWriter == null)
            return;

        for(JunctionData junctionData : bucket.junctions())
        {
            junctionData.JunctionGroups.forEach(x -> mBamWriter.writeRecords(x.reads()));
            mBamWriter.writeRecords(junctionData.SupportingReads);
        }
    }
}
