package com.hartwig.hmftools.svprep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
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
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.svprep.reads.JunctionData;
import com.hartwig.hmftools.svprep.reads.ReadFilterType;
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.RemoteJunction;

public class ResultsWriter
{
    private final SvConfig mConfig;

    private final BufferedWriter mReadWriter;
    private final BufferedWriter mBedWriter;
    private final BufferedWriter mJunctionWriter;
    private final BamWriter mBamWriter;

    public ResultsWriter(final SvConfig config)
    {
        mConfig = config;

        if(mConfig.OutputDir == null)
        {
            mBedWriter = null;
            mReadWriter = null;
            mJunctionWriter = null;
            mBamWriter = null;
            return;
        }

        mJunctionWriter = initialiseJunctionWriter();
        mBedWriter = initialiseBedWriter();
        mReadWriter = initialiseReadWriter();
        mBamWriter = new BamWriter(config);
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mBedWriter);
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

            writer.write("ReadId,GroupCount,GroupStatus,HasExternal,ReadType,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",FragLength,MateChr,MatePosStart,MapQual,MultiMapped,SuppData,Flags");
            writer.write(",FirstInPair,ReadReversed,Proper,Duplicate,Unmapped,MateUnmapped,Secondary,Supplementary,JunctionPositions");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeReadData(final List<ReadGroup> readGroups)
    {
        if(mReadWriter == null)
            return;

        try
        {
            for(ReadGroup readGroup : readGroups)
            {
                String junctionPosStr = "";

                if(readGroup.junctionPositions() != null)
                {
                    StringJoiner sjPos = new StringJoiner(ITEM_DELIM);
                    readGroup.junctionPositions().forEach(x -> sjPos.add(String.valueOf(x)));
                    junctionPosStr = sjPos.toString();
                }

                for(ReadRecord read : readGroup.reads())
                {
                    mReadWriter.write(format("%s,%d,%s,%s",
                            read.id(), readGroup.size(), readGroup.groupStatus(), readGroup.spansPartitions()));

                    mReadWriter.write(format(",%s,%s,%d,%d,%s",
                            read.readType(), read.Chromosome, read.start(), read.end(), read.cigar().toString()));

                    SupplementaryReadData suppData = read.supplementaryAlignment();

                    mReadWriter.write(format(",%d,%s,%d,%d,%s,%s,%d",
                            read.fragmentInsertSize(), read.MateChromosome, read.MatePosStart, read.MapQuality,
                            read.isMultiMapped(), suppData != null ? suppData.asCsv() : "N/A", read.flags()));

                    mReadWriter.write(format(",%s,%s,%s,%s,%s,%s,%s,%s",
                            read.isFirstOfPair(), read.isReadReversed(),
                            read.hasFlag(PROPER_PAIR), read.isDuplicate(), read.hasFlag(READ_UNMAPPED), read.hasFlag(MATE_UNMAPPED),
                            read.hasFlag(SECONDARY_ALIGNMENT), read.hasFlag(SUPPLEMENTARY_ALIGNMENT)));

                    mReadWriter.write(format(",%s", junctionPosStr));

                    mReadWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }

    private BufferedWriter initialiseJunctionWriter()
    {
        if(!mConfig.WriteTypes.contains(JUNCTIONS))
            return null;

        try
        {
            String filename = mConfig.formFilename(JUNCTIONS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Chromosome,Position,Orientation,Fragments,SupportingReads,LowMapQualFrags,Hotspot,InitialReadId");
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

    public synchronized void writeBamRecords(final List<ReadGroup> readGroups)
    {
        readGroups.forEach(x -> writeBamRecords(x));
    }

    public synchronized void writeBamRecords(final ReadGroup readGroup)
    {
        if(mBamWriter == null)
            return;

        mBamWriter.writeRecords(readGroup.reads());
    }

    private BufferedWriter initialiseBedWriter()
    {
        if(!mConfig.WriteTypes.contains(SV_BED))
            return null;

        return null;
    }
}
