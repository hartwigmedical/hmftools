package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.BAM_RECORD_SAMPLE_ID_TAG;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_EXACT_SUPPORT_FRAGS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_EXTRA_INFO;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_HOTSPOT_JUNCTION;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_INDEL_JUNCTION;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_JUNCTION_FRAGS;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FLD_OTHER_SUPPORT_FRAGS;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.JUNCTIONS;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.READS;

import static htsjdk.samtools.SAMFlag.MATE_REVERSE_STRAND;
import static htsjdk.samtools.SAMFlag.READ_UNMAPPED;
import static htsjdk.samtools.SAMFlag.SUPPLEMENTARY_ALIGNMENT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;
import com.hartwig.hmftools.esvee.prep.types.RemoteJunction;

public class ResultsWriter
{
    private final PrepConfig mConfig;

    private final BufferedWriter mReadWriter;
    private final BufferedWriter mJunctionWriter;
    private final BamWriter mBamWriter;

    public ResultsWriter(final PrepConfig config)
    {
        mConfig = config;

        if(mConfig.OutputDir == null)
        {
            mReadWriter = null;
            mJunctionWriter = null;
            mBamWriter = null;
            return;
        }

        mJunctionWriter = initialiseJunctionWriter();
        mReadWriter = initialiseReadWriter();
        mBamWriter = new BamWriter(config);
    }

    public void close()
    {
        closeBufferedWriter(mReadWriter);
        closeBufferedWriter(mJunctionWriter);
        mBamWriter.close();
    }

    public long writtenCount() { return mBamWriter != null ? mBamWriter.writtenCount() : 0; }

    public synchronized void writeReadGroup(final List<ReadGroup> readGroups)
    {
        for(ReadGroup readGroup : readGroups)
        {
            if(filterReadGroup(readGroup))
                continue;

            if(mReadWriter != null)
            {
                String junctionPosStr = readGroup.junctionPositionsStr();

                for(PrepRead read : readGroup.reads())
                {
                    if(read.written())
                        continue;

                    writeReadData(
                            read, readGroup.size(), readGroup.expectedReadCount(), readGroup.groupStatus(), readGroup.spansPartitions(),
                            junctionPosStr);
                }
            }

            writeBamRecords(readGroup);

            readGroup.reads().forEach(x -> x.setWritten());
        }
    }

    private static boolean filterReadGroup(final ReadGroup readGroup)
    {
        if(readGroup.conditionalOnRemoteReads() && !readGroup.hasRemoteJunctionReads())
            return true;

        return false;
    }

    private BufferedWriter initialiseJunctionWriter()
    {
        if(!mConfig.WriteTypes.contains(JUNCTIONS))
            return null;

        try
        {
            String filename = mConfig.formFilename(JUNCTIONS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_POSITION);
            sj.add(FLD_ORIENTATION);
            sj.add(FLD_JUNCTION_FRAGS);
            sj.add(FLD_EXACT_SUPPORT_FRAGS);
            sj.add(FLD_OTHER_SUPPORT_FRAGS);
            sj.add("LowMapQualFrags");
            sj.add("MaxQual");
            sj.add(FLD_EXTRA_INFO);
            sj.add(FLD_INDEL_JUNCTION);
            sj.add(FLD_HOTSPOT_JUNCTION);
            sj.add("InitialReadId");

            if(mConfig.TrackRemotes)
            {
                sj.add("RemoteJunctionCount");
                sj.add("RemoteJunctions");
            }

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to create junction writer: {}", e.toString());
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
                int extraInfoValues = 0;

                int junctionFrags = junctionData.junctionFragmentCount();
                int exactSupportFrags = junctionData.ExactSupportGroups.size();
                int otherSupportFrags = junctionData.SupportingGroups.size();

                if(!junctionData.discordantGroup())
                {
                    boolean expectLeftClipped = junctionData.Orient.isReverse();

                    for(PrepRead read : junctionData.ReadTypeReads.get(ReadType.JUNCTION))
                    {
                        // check the read supports this junction (it can also support another junction)
                        boolean supportsJunction =
                                (expectLeftClipped && read.start() == junctionData.Position && read.isLeftClipped())
                                        || (!expectLeftClipped && read.end() == junctionData.Position && read.isRightClipped());

                        if(!supportsJunction)
                            continue;

                        if(ReadFilterType.isSet(read.filters(), ReadFilterType.MIN_MAP_QUAL))
                            ++lowMapQualFrags;

                        maxMapQual = max(maxMapQual, read.mapQuality());

                        if(!junctionData.internalIndel())
                        {
                            int scLength = expectLeftClipped ? read.leftClipLength() : read.rightClipLength();

                            if(scLength > extraInfoValues)
                            {
                                extraInfoValues = scLength;
                            }
                        }
                    }

                    for(PrepRead read : junctionData.ReadTypeReads.get(ReadType.EXACT_SUPPORT))
                    {
                        maxMapQual = max(maxMapQual, read.mapQuality());

                        if(ReadFilterType.isSet(read.filters(), ReadFilterType.MIN_MAP_QUAL))
                            ++lowMapQualFrags;
                    }
                }
                else
                {
                    // replace soft-clip length with max remote location reads
                    extraInfoValues = !junctionData.RemoteJunctions.isEmpty() ? junctionData.RemoteJunctions.get(0).Fragments : 0;
                }

                mJunctionWriter.write(String.format("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                        chromosome, junctionData.Position, junctionData.Orient.asByte(), junctionFrags,
                        exactSupportFrags, otherSupportFrags, lowMapQualFrags, maxMapQual));

                mJunctionWriter.write(String.format("\t%d\t%s\t%s\t%s",
                        extraInfoValues, junctionData.internalIndel(), junctionData.hotspot(),
                        junctionData.topJunctionRead() != null ? junctionData.topJunctionRead().id() : "EXISTING"));

                if(mConfig.TrackRemotes)
                {
                    // RemoteChromosome:RemotePosition:RemoteOrientation;Fragments then separated by ';'
                    String remoteJunctionsStr = "";

                    if(!junctionData.RemoteJunctions.isEmpty())
                    {
                        Collections.sort(junctionData.RemoteJunctions, Comparator.comparingInt(x -> -x.Fragments));

                        StringJoiner sj = new StringJoiner(ITEM_DELIM);

                        for(int i = 0; i < min(junctionData.RemoteJunctions.size(), 5); ++i)
                        {
                            RemoteJunction remoteJunction = junctionData.RemoteJunctions.get(i);
                            sj.add(String.format("%s:%d:%d:%d",
                                    remoteJunction.Chromosome, remoteJunction.Position, remoteJunction.Orient.asByte(), remoteJunction.Fragments));
                        }
                        remoteJunctionsStr = sj.toString();
                    }

                    mJunctionWriter.write(String.format("\t%d\t%s", junctionData.RemoteJunctions.size(), remoteJunctionsStr));
                }

                mJunctionWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write junction data: {}", e.toString());
        }
    }

    private void writeBamRecords(final ReadGroup readGroup)
    {
        if(mBamWriter == null)
            return;

        // note additional filters for a read to be written to the BAM
        // - excessive low qual soft-clip bases
        // - above the poly-G(C) threshold
        // - cannot be a group of only supplementaries (in case the group is an unmarked duplicate)
        for(PrepRead read : readGroup.reads())
        {
            if(filterBamRecord(read))
                continue;

            if(read.written())
                continue;

            mBamWriter.writeRecord(read.record());
        }
    }

    private boolean filterBamRecord(final PrepRead read)
    {
        if(ReadFilterType.isSet(read.filters(), ReadFilterType.POLY_G_SC))
            return true;

        return false;
    }

    private BufferedWriter initialiseReadWriter()
    {
        if(!mConfig.WriteTypes.contains(READS))
            return null;

        try
        {
            String filename = mConfig.formFilename(READS);
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("ReadId").add("SampleId").add("GroupCount").add("ExpectedCount").add("GroupStatus").add("HasExternal").add("ReadType");
            sj.add("Chromosome").add("PosStart").add("PosEnd").add("Cigar");
            sj.add("FragLength").add("MateChr").add("MatePosStart").add("MapQual").add("SuppData").add("Flags").add("Filters");
            sj.add("FirstInPair").add("ReadReversed").add("MateReversed").add("Unmapped").add("MateUnmapped").add("Supplementary");
            sj.add("JunctionPositions");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to create read writer: {}", e.toString());
        }

        return null;
    }

    private void writeReadData(
            final PrepRead read, int readCount, int expectedReadCount, final ReadGroupStatus status, boolean spansPartitions,
            final String junctionPositions)
    {
        if(mReadWriter == null)
            return;

        try
        {
            String sampleId = read.record().getStringAttribute(BAM_RECORD_SAMPLE_ID_TAG);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(read.id());
            sj.add(sampleId);
            sj.add(String.valueOf(readCount));
            sj.add(String.valueOf(expectedReadCount));
            sj.add(String.valueOf(status));
            sj.add(String.valueOf(spansPartitions));
            sj.add(String.valueOf(read.readType()));
            sj.add(read.Chromosome);
            sj.add(String.valueOf(read.start()));
            sj.add(String.valueOf(read.end()));
            sj.add(String.valueOf(read.cigar()));
            sj.add(String.valueOf(read.fragmentInsertSize()));
            sj.add(read.MateChromosome);

            sj.add(String.valueOf(read.MatePosStart));
            sj.add(String.valueOf(read.mapQuality()));

            SupplementaryReadData suppData = read.supplementaryAlignment();
            sj.add(suppData != null ? suppData.asDelimStr() : "N/A");

            sj.add(String.valueOf(read.flags()));
            sj.add(String.valueOf(read.filters()));
            sj.add(String.valueOf(read.isFirstOfPair()));
            sj.add(String.valueOf(read.isReadReversed()));
            sj.add(String.valueOf(read.hasFlag(MATE_REVERSE_STRAND)));
            sj.add(String.valueOf(read.hasFlag(READ_UNMAPPED)));
            sj.add(String.valueOf(read.hasMate() && read.isMateUnmapped()));
            sj.add(String.valueOf(read.hasFlag(SUPPLEMENTARY_ALIGNMENT)));
            sj.add(String.valueOf(junctionPositions));

            mReadWriter.write(sj.toString());
            mReadWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error(" failed to write read data: {}", e.toString());
        }
    }
}
