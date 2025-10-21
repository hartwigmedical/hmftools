package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundDuplicateCount;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundFragmentLength;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundGcPercent;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.TARGET_REGION_PROXIMITY;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_INFO_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.sv.SvUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionReader extends Thread
{
    private final MetricsConfig mConfig;

    private final Queue<Partition> mPartitions;
    private final int mPartitionCount;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
//    private final BufferedWriter mReadDataWriter;

    private Partition mCurrentPartition;

//    private final FragmentGcMap mAllFragmentGcMap; // will contribute to whole BAM metrics
//    private final FragmentGcMap mNonTargetedFragmentGcMap; // for reads outside targeted reginos
//    private final FragmentGcMap mTargetedFragmentGcMap;

    private final Map<String,SAMRecord> mReadGroupMap;

//    List<SAMRecord> samRecordsOfInterest = new ArrayList<>();

    public PartitionReader(final MetricsConfig config, final Queue<Partition> partitions, final BufferedWriter readDataWriter)
    {
        mConfig = config;

        mPartitions = partitions;
        mPartitionCount = partitions.size();
//        mReadDataWriter = readDataWriter;
        mCurrentPartition = null;

        mSamReader = SamReaderFactory.makeDefault().open(new File(mConfig.BamFile));

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false, false, false);
        mBamSlicer.setKeepUnmapped();

        mReadGroupMap = Maps.newHashMap();
//        mAllFragmentGcMap = new FragmentGcMap();
//        mNonTargetedFragmentGcMap = new FragmentGcMap();
//        mTargetedFragmentGcMap = new FragmentGcMap();
    }

    @Override
    public void run()
    {
        if(mPartitions.isEmpty())
            return;

        while(true)
        {
            try
            {
                mCurrentPartition = mPartitions.remove();

                int processedCount = mPartitionCount - mPartitions.size();

                if((processedCount % 100) == 0)
                {
                    CB_LOGGER.info("processed {} partitions", processedCount);
                }

                processPartition();
            }
            catch(NoSuchElementException e)
            {
                CB_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }

//    public FragmentGcMap allFragmentGcMap() { return mAllFragmentGcMap; }
//    public FragmentGcMap nonTargetedFragmentGcMap() { return mNonTargetedFragmentGcMap; }
//    public FragmentGcMap targetedFragmentGcMap() { return mTargetedFragmentGcMap; }

    private void processPartition()
    {
        CB_LOGGER.trace("processing region({})", mCurrentPartition);
        mBamSlicer.slice(mSamReader, mCurrentPartition, this::processSamRecord);
        CB_LOGGER.trace("completed region({})", mCurrentPartition);
        postSliceProcess();
    }

    @VisibleForTesting
    protected void postSliceProcess()
    {
        // process overlapping groups
        for(SAMRecord read : mReadGroupMap.values())
        {
            processSingleRecord(read);
        }

        mReadGroupMap.clear();

//        if (!samRecordsOfInterest.isEmpty())
//        {
//            System.out.println("Processing " + samRecordsOfInterest.size() + " SAM records");
//        }
    }

    private void processSamRecord(final SAMRecord read)
    {
        processSamRecord(read, true);
    }

//    private void captureRecordIfIntersectsTargetRegion(SAMRecord read)
//    {
//        int fragBoundsStart = read.getAlignmentStart() - TARGET_REGION_PROXIMITY;
//        int fragBoundsEnd = read.getAlignmentEnd()+ TARGET_REGION_PROXIMITY;
//
//        for(TargetRegionData targetRegion : mCurrentPartition.TargetRegions)
//        {
//            if(positionsOverlap(fragBoundsStart, fragBoundsEnd, targetRegion.start(), targetRegion.end()))
//            {
//                samRecordsOfInterest.add(read);
//            }
//        }
////        samRecordsOfInterest.add(read);
//
//    }

    private void processSamRecord(final SAMRecord read, boolean removeFragments)
    {
//        captureRecordIfIntersectsTargetRegion(read);
        int readStart = read.getAlignmentStart();

        if(!mCurrentPartition.containsPosition(readStart))
            return;

        if(SvUtils.isDiscordant(read))
        {
            processSingleRecord(read);
            return;
        }

        SAMRecord mateRead = mReadGroupMap.get(read.getReadName());

        if(mateRead == null)
        {
            mReadGroupMap.put(read.getReadName(), read);
            return;
        }

        processFragment(read, mateRead);

        if(removeFragments)
            mReadGroupMap.remove(read.getReadName());
    }

    private void processSingleRecord(final SAMRecord read)
    {
        int duplicateCount = getDuplicateReadCount(read);

        String alignedBases = getAlignedReadBases(read, 0);
        double gcContent = GcCalcs.calcGcPercent(alignedBases);

        int readPosStart = read.getAlignmentStart();
        int readPosEnd = read.getAlignmentEnd();

        if(read.getReadNegativeStrandFlag())
            readPosEnd += CigarUtils.rightSoftClipLength(read);
        else
            readPosStart -= CigarUtils.leftSoftClipLength(read);

        addFragmentData(
                read.getReadName(), 1, readPosStart, readPosEnd, 0, alignedBases.length(), gcContent, duplicateCount);
    }

    private void processFragment(final SAMRecord read, final SAMRecord mate)
    {
        if(read.getMappingQuality() < MIN_MAPPING_QUALITY || mate.getMappingQuality() < MIN_MAPPING_QUALITY)
            return;

        int duplicateCount = getDuplicateReadCount(read);

        int readPosStart = read.getAlignmentStart();
        int readPosEnd = read.getAlignmentEnd();

        if(read.getReadNegativeStrandFlag())
            readPosEnd += CigarUtils.rightSoftClipLength(read);
        else
            readPosStart -= CigarUtils.leftSoftClipLength(read);

        int matePosStart = mate.getAlignmentStart();
        int matePosEnd = mate.getAlignmentEnd();

        if(mate.getReadNegativeStrandFlag())
            matePosEnd += CigarUtils.rightSoftClipLength(mate);
        else
            matePosStart -= CigarUtils.leftSoftClipLength(mate);

        int fragPosStart = min(readPosStart, matePosStart);
        int fragPosEnd = max(readPosEnd, matePosEnd);

        double gcContent;
        String alignedBases;

        if(positionsOverlap(readPosStart, readPosEnd, matePosStart, matePosEnd))
        {
            SAMRecord lowerRead = readPosStart <= matePosStart ? read : mate;
            SAMRecord upperRead = read == lowerRead ? mate : read;

            alignedBases = getAlignedReadBases(lowerRead, 0);

            int minReadEnd = min(readPosEnd, matePosEnd);
            int maxReadEnd = max(readPosEnd, matePosEnd);

            if(maxReadEnd > minReadEnd)
            {
                alignedBases += getAlignedReadBases(upperRead, minReadEnd + 1);
            }

            gcContent = GcCalcs.calcGcPercent(alignedBases);
        }
        else
        {
            String fragmentBases = getAlignedReadBases(read, 0);
            fragmentBases += getAlignedReadBases(mate, 0);
            alignedBases = fragmentBases;

            int gapPosStart = readPosEnd < matePosStart ? readPosEnd : matePosEnd + 1;
            int gapPosEnd = readPosEnd < matePosStart ? matePosStart - 1 : readPosStart - 1;

//            String gapRefBases = mConfig.RefGenome.getBaseString(read.getReferenceName(), gapPosStart, gapPosEnd);
//            fragmentBases += gapRefBases;

            gcContent = GcCalcs.calcGcPercent(fragmentBases);
        }

        int insertSize = abs(read.getInferredInsertSize());

        addFragmentData(read.getReadName(), 2, fragPosStart, fragPosEnd, insertSize, alignedBases.length(), gcContent, duplicateCount);
    }

    private static int getDuplicateReadCount(final SAMRecord read)
    {
        Object consensusReadAttr = read.getAttribute(CONSENSUS_READ_ATTRIBUTE);

        if(consensusReadAttr == null)
            return 1;

        if(consensusReadAttr instanceof String)
            return Integer.parseInt(String.valueOf(consensusReadAttr).split(CONSENSUS_INFO_DELIM, 2)[0]);

        return Integer.parseInt(consensusReadAttr.toString());
    }

    private void addFragmentData(
            final String readId, int readCount, int fragmentPosStart, int fragmentPosEnd, int rawFragmentLength, int baseCount,
            double rawGcPercent, int duplicateCount)
    {
        mCurrentPartition.recordFragment(fragmentPosStart, rawFragmentLength);

//        double gcPercent = roundGcPercent(rawGcPercent, mConfig.GcPercentUnits);
//        int fragmentLength = fragmentPosEnd - fragmentPosStart + 1;
//        int roundedFragmentLength = roundFragmentLength(rawFragmentLength, mConfig.FragmentLengthUnits);
//
//        int dupCountRounded = roundDuplicateCount(duplicateCount);

//        mAllFragmentGcMap.add(roundedFragmentLength, gcPercent, dupCountRounded);

//        boolean matchesRegion = false;
//
//        int fragBoundsStart = fragmentPosStart - TARGET_REGION_PROXIMITY;
//        int fragBoundsEnd = fragmentPosEnd + TARGET_REGION_PROXIMITY;
//
//        for(TargetRegionData targetRegion : mCurrentPartition.TargetRegions)
//        {
//            if(positionsOverlap(fragBoundsStart, fragBoundsEnd, targetRegion.start(), targetRegion.end()))
//            {
//                matchesRegion = true;
//
//                if(mConfig.CaptureRegionCounts)
//                    targetRegion.FragmentGcCounts.add(roundedFragmentLength, gcPercent, dupCountRounded);
//                else
//                    break;
//            }
//        }

//        if(matchesRegion)
//            mTargetedFragmentGcMap.add(roundedFragmentLength, gcPercent, dupCountRounded);
//        else
//            mNonTargetedFragmentGcMap.add(roundedFragmentLength, gcPercent, dupCountRounded);
//
//        if(mReadDataWriter != null)
//        {
//            PartitionReader.writeReadData(
//                    mReadDataWriter, readId, readCount, mCurrentPartition.Chromosome, fragmentPosStart, fragmentPosEnd,
//                    rawFragmentLength, fragmentLength, baseCount, rawGcPercent);
//        }
    }

    private static String getAlignedReadBases(final SAMRecord read, final int minReadStartPos)
    {
        String readBasesStr = read.getReadString();

        if(read.getCigar().getCigarElements().size() == 1)
        {
            if(minReadStartPos > 0)
            {
                int readStartOffset = minReadStartPos - read.getAlignmentStart();
                return readBasesStr.substring(readStartOffset);
            }
            else
            {
                return readBasesStr;
            }
        }

        boolean useLeftClip = minReadStartPos == 0 && !read.getReadNegativeStrandFlag();
        boolean useRightClip = read.getReadNegativeStrandFlag();

        int readIndex = 0;
        int position = read.getAlignmentStart();

        String alignedReadBases = "";

        for(int i = 0; i < read.getCigar().getCigarElements().size(); ++i)
        {
            CigarElement element = read.getCigar().getCigarElements().get(i);

            switch(element.getOperator())
            {
                case S:

                    if((i == 0 && useLeftClip) || (i > 0 && useRightClip))
                    {
                        int readIndexStart = readIndex;
                        int readIndexEnd = readIndexStart + element.getLength() - 1;
                        alignedReadBases += readBasesStr.substring(readIndexStart, readIndexEnd);
                    }

                    readIndex += element.getLength();
                    break;

                case D:
                case N:
                    position += element.getLength();
                    break;

                case H:
                    break;

                case I:
                    readIndex += element.getLength();
                    break;

                case M:

                    int alignedEnd = position + element.getLength() - 1;

                    if(minReadStartPos <= alignedEnd)
                    {
                        // eg pos = 5, min read pos = 10, aligned end == 12
                        // then offset = 5, read bases are current read index + 5, to 5 + 8 - 5
                        int offsetStart = max(minReadStartPos - position, 0);
                        int readIndexStart = readIndex + offsetStart;
                        int readIndexEnd = readIndexStart + element.getLength() - offsetStart;
                        alignedReadBases += readBasesStr.substring(readIndexStart, readIndexEnd);
                    }

                    position += element.getLength();
                    readIndex += element.getLength();
                    break;

                default:
                    break;
            }
        }

        return alignedReadBases;
    }

//    public static BufferedWriter initialiseReadWriter(final MetricsConfig config)
//    {
//        try
//        {
//            // write summary metrics
//            String filename = config.OutputDir + config.SampleId + ".read_data";
//
//            if(config.OutputId != null)
//                filename += "." + config.OutputId;
//
//            filename += TSV_EXTENSION;
//
//            BufferedWriter writer = createBufferedWriter(filename);
//
//            StringJoiner sj = new StringJoiner(TSV_DELIM);
//            sj.add("ReadId").add("ReadCount").add("Chromosome").add("FragmentStart").add("FragmentEnd").add("RawInsertSize");
//            sj.add("FragmentLength").add("BaseCount").add("GcRatio");
//            writer.write(sj.toString());
//            writer.newLine();
//            return writer;
//
//        }
//        catch(IOException e)
//        {
//            CB_LOGGER.error("failed to initialise read data writer: {}", e.toString());
//            return null;
//        }
//    }

    public synchronized static void writeReadData(
            final BufferedWriter writer, final String readId, int readCount, final String chromosome, int fragmentStart, int fragmentEnd,
            int rawInsertSize, int fragmentLength, int baseCount, double gcRatio)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(readId);
            sj.add(String.valueOf(readCount));
            sj.add(chromosome);
            sj.add(String.valueOf(fragmentStart));
            sj.add(String.valueOf(fragmentEnd));
            sj.add(String.valueOf(rawInsertSize));
            sj.add(String.valueOf(fragmentLength));
            sj.add(String.valueOf(baseCount));
            sj.add(format("%.3f", gcRatio));
            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to write read data: {}", e.toString());
        }
    }

}
