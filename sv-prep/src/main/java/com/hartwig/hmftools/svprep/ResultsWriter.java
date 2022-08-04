package com.hartwig.hmftools.svprep;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.JUNCTIONS;
import static com.hartwig.hmftools.svprep.WriteType.READS;
import static com.hartwig.hmftools.svprep.WriteType.SV_BED;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.MIN_MAP_QUAL;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.POLY_G_SC;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.SOFT_CLIP_LOW_BASE_QUAL;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.hasPolyATSoftClip;

import static htsjdk.samtools.SAMFlag.DUPLICATE_READ;
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
import com.hartwig.hmftools.svprep.reads.ReadGroup;
import com.hartwig.hmftools.svprep.reads.ReadGroupStatus;
import com.hartwig.hmftools.svprep.reads.ReadRecord;
import com.hartwig.hmftools.svprep.reads.ReadType;
import com.hartwig.hmftools.svprep.reads.RemoteJunction;

import htsjdk.samtools.SAMRecord;

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

            writer.write("ReadId,GroupCount,ExpectedCount,GroupStatus,HasExternal,ReadType,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",FragLength,MateChr,MatePosStart,MapQual,SuppData,Flags,Filters");
            writer.write(",FirstInPair,ReadReversed,Proper,Unmapped,MateUnmapped,Supplementary,Duplicate,JunctionPositions");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    public synchronized void writeReadGroup(final List<ReadGroup> readGroups)
    {
        for(ReadGroup readGroup : readGroups)
        {
            if(readGroup.onlySupplementaries() && !readGroup.hasRemoteNonSupplementaries())
                continue;

            String junctionPosStr = "";

            if(readGroup.junctionPositions() != null)
            {
                StringJoiner sjPos = new StringJoiner(ITEM_DELIM);
                readGroup.junctionPositions().forEach(x -> sjPos.add(String.valueOf(x)));
                junctionPosStr = sjPos.toString();
            }

            for(ReadRecord read : readGroup.reads())
            {
                if(read.written())
                    continue;

                writeReadData(
                        read, readGroup.size(), readGroup.expectedReadCount(), readGroup.groupStatus(), readGroup.spansPartitions(),
                        junctionPosStr);
            }
        }
    }

    private void writeReadData(
            final ReadRecord read, int readCount, int expectedReadCount, final ReadGroupStatus status, boolean spansPartitions,
            final String junctionPositions)
    {
        if(mReadWriter == null || read.written())
            return;

        try
        {
            mReadWriter.write(format("%s,%d,%d,%s,%s", read.id(), readCount, expectedReadCount, status, spansPartitions));

            mReadWriter.write(format(",%s,%s,%d,%d,%s",
                    read.readType(), read.Chromosome, read.start(), read.end(), read.cigar().toString()));

            SupplementaryReadData suppData = read.supplementaryAlignment();

            mReadWriter.write(format(",%d,%s,%d,%d,%s,%d,%d",
                    read.fragmentInsertSize(), read.MateChromosome, read.MatePosStart, read.mapQuality(),
                    suppData != null ? suppData.asCsv() : "N/A", read.flags(), read.filters()));

            mReadWriter.write(format(",%s,%s,%s,%s,%s,%s,%s",
                    read.isFirstOfPair(), read.isReadReversed(), read.hasFlag(PROPER_PAIR), read.hasFlag(READ_UNMAPPED),
                    read.isMateUnmapped(), read.hasFlag(SUPPLEMENTARY_ALIGNMENT), read.hasFlag(DUPLICATE_READ)));

            mReadWriter.write(format(",%s", junctionPositions));

            mReadWriter.newLine();
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

            writer.write("Chromosome,Position,Orientation,JunctionFrags,SupportFrags,DiscordantFrags,LowMapQualFrags,MaxQual");
            writer.write(",MaxSoftClip,BaseDepth,HasPolyAT,Indel,Hotspot,SoftClipBases,InitialReadId,RemoteJunctionCount,RemoteJunctions");
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
                int maxMapQual = 0;
                int lowMapQualFrags = 0;
                int maxSoftClip = 0;
                String softClipBases = "";
                boolean hasPloyAT = false;
                boolean expectLeftClipped = junctionData.Orientation == NEG_ORIENT;

                for(ReadRecord read : junctionData.ReadTypeReads.get(ReadType.JUNCTION))
                {
                    // check the read supports this junction (it can only support another junction)
                    boolean supportsJunction =
                            (expectLeftClipped && read.start() == junctionData.Position && read.cigar().isLeftClipped())
                            || (!expectLeftClipped && read.end() == junctionData.Position && read.cigar().isRightClipped());

                    if(!supportsJunction)
                        continue;

                    if(ReadFilterType.isSet(read.filters(), MIN_MAP_QUAL))
                        ++lowMapQualFrags;

                    maxMapQual = max(maxMapQual, read.mapQuality());

                    if(!junctionData.internalIndel())
                    {
                        if(!hasPloyAT)
                            hasPloyAT = hasPolyATSoftClip(read, expectLeftClipped);

                        int scLength = expectLeftClipped ?
                                read.cigar().getFirstCigarElement().getLength() : read.cigar().getLastCigarElement().getLength();

                        if(scLength > maxSoftClip)
                        {
                            maxSoftClip = scLength;
                            softClipBases = ReadRecord.getSoftClippedBases(read.record(), expectLeftClipped);
                        }
                    }
                }

                int exactSupportFrags = junctionData.ExactSupportGroups.size();
                int discordantFrags = junctionData.SupportingGroups.size();

                for(ReadRecord read : junctionData.ReadTypeReads.get(ReadType.EXACT_SUPPORT))
                {
                    ++exactSupportFrags;
                    maxMapQual = max(maxMapQual, read.mapQuality());

                    if(ReadFilterType.isSet(read.filters(), MIN_MAP_QUAL))
                        ++lowMapQualFrags;
                }

                mJunctionWriter.write(format("%s,%d,%d,%d,%d,%d,%d,%d",
                        chromosome, junctionData.Position, junctionData.Orientation, junctionData.junctionFragmentCount(),
                        exactSupportFrags, discordantFrags, lowMapQualFrags, maxMapQual));

                mJunctionWriter.write(format(",%d,%d,%s,%s,%s,%s,%s",
                        maxSoftClip, junctionData.depth(), hasPloyAT, junctionData.internalIndel(), junctionData.hotspot(),
                        softClipBases, junctionData.InitialRead != null ? junctionData.InitialRead.id() : "EXISTING"));

                // RemoteChromosome:RemotePosition:RemoteOrientation;Fragments then separated by ';'
                String remoteJunctionsStr = "";

                if(!junctionData.RemoteJunctions.isEmpty())
                {
                    Collections.sort(junctionData.RemoteJunctions, new RemoteJunction.RemoteJunctionSorter());

                    StringJoiner sj = new StringJoiner(ITEM_DELIM);

                    for(int i = 0; i < min(junctionData.RemoteJunctions.size(), 10); ++i)
                    {
                        RemoteJunction remoteJunction = junctionData.RemoteJunctions.get(i);
                        sj.add(format("%s:%d:%d;%d",
                                remoteJunction.Chromosome, remoteJunction.Position, remoteJunction.Orientation, remoteJunction.Fragments));
                        // junctionData.RemoteJunctions.forEach(x -> sj.add(format("%s:%d:%d", x.Chromosome, x.Position, x.Orientation)));
                    }
                    remoteJunctionsStr = sj.toString();
                }

                mJunctionWriter.write(format(",%d,%s", junctionData.RemoteJunctions.size(), remoteJunctionsStr));
                mJunctionWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write junction data: {}", e.toString());
        }
    }

    public synchronized void writeBamRecords(final List<ReadGroup> readGroups)
    {
        if(mBamWriter == null)
            return;

        // note additional filters for a read to be written to the BAM
        // - excessive low qual soft-clip bases
        // - above the poly-G(C) threshold
        // - cannot be a group of only supplementaries (in case the group is an unmarked duplicate)

        for(ReadGroup readGroup : readGroups)
        {
            if(readGroup.onlySupplementaries() && !readGroup.hasRemoteNonSupplementaries())
                continue;

            readGroup.reads().stream().filter(x -> !filterBamRecord(x)).forEach(x -> mBamWriter.writeRecord(x.record()));
        }
    }

    private boolean filterBamRecord(final ReadRecord read)
    {
        if(read.written())
            return true;

        if(ReadFilterType.isSet(read.filters(), POLY_G_SC))
            return true;

        if(ReadFilterType.isSet(read.filters(), SOFT_CLIP_LOW_BASE_QUAL))
        {
            if(read.readType() != ReadType.JUNCTION && read.readType() != ReadType.EXACT_SUPPORT)
                return true;
        }

        return false;
    }

    private BufferedWriter initialiseBedWriter()
    {
        if(!mConfig.WriteTypes.contains(SV_BED))
            return null;

        return null;
    }
}
