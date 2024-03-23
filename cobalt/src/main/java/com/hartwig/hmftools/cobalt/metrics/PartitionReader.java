package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundDuplicateCount;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundFragmentLength;
import static com.hartwig.hmftools.cobalt.metrics.FragmentGcCounts.roundGcPercent;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.TARGET_REGION_PROXIMITY;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_INFO_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
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

    private Partition mCurrentPartition;

    private final FragmentGcMap mAllFragmentGcMap; // will contribute to whole BAM metrics
    private final FragmentGcMap mNonTargetedFragmentGcMap; // for reads outside targeted reginos
    private final FragmentGcMap mTargetedFragmentGcMap;

    private final Map<String,SAMRecord> mReadGroupMap;

    public PartitionReader(final MetricsConfig config, final Queue<Partition> partitions)
    {
        mConfig = config;

        mPartitions = partitions;
        mPartitionCount = partitions.size();
        mCurrentPartition = null;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false, false, false);
        mBamSlicer.setKeepUnmapped();

        mReadGroupMap = Maps.newHashMap();
        mAllFragmentGcMap = new FragmentGcMap();
        mNonTargetedFragmentGcMap = new FragmentGcMap();
        mTargetedFragmentGcMap = new FragmentGcMap();
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

    public FragmentGcMap allFragmentGcMap() { return mAllFragmentGcMap; }
    public FragmentGcMap nonTargetedFragmentGcMap() { return mNonTargetedFragmentGcMap; }
    public FragmentGcMap targetedFragmentGcMap() { return mTargetedFragmentGcMap; }

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
            processSamRecord(read);
        }

        mReadGroupMap.clear();
    }

    private void processSamRecord(final SAMRecord read)
    {
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
    }

    private void processSingleRecord(final SAMRecord read)
    {
        int duplicateCount = getDuplicateReadCount(read);
        double gcContent = getReadGcContent(read);

        addFragmentData(read.getAlignmentStart(), read.getAlignmentEnd(), duplicateCount, gcContent, 0);
    }

    private void processFragment(final SAMRecord read, final SAMRecord mate)
    {
        int duplicateCount = getDuplicateReadCount(read);

        int fragmentLength = abs(read.getInferredInsertSize());

        // CHECK: should soft-clipped bases be included or ignored?
        int readPosStart = read.getAlignmentStart();
        int readPosEnd = read.getAlignmentEnd();

        int matePosStart = mate.getAlignmentStart();
        int matePosEnd = mate.getAlignmentEnd();

        int fragPosStart = min(readPosStart, matePosStart);
        int fragPosEnd = max(readPosEnd, matePosEnd);

        double gcContent;

        if(positionsOverlap(readPosStart, readPosEnd, matePosStart, matePosEnd))
        {
            // CHECK: how accurate to make this for overlapping fragments? for now just take the lower read's bases
            // and append with GC content of upper non-overlapping bases
            SAMRecord lowerRead = read.getAlignmentStart() <= mate.getAlignmentStart() ? read : mate;
            SAMRecord upperRead = read == lowerRead ? mate : read;

            String alignedBases = getAlignedReadBases(lowerRead, 0);

            if(upperRead.getAlignmentEnd() > lowerRead.getAlignmentEnd())
            {
                alignedBases += getAlignedReadBases(upperRead, lowerRead.getAlignmentEnd());
            }

            gcContent = GcCalcs.calcGcPercent(alignedBases);
        }
        else
        {
            String fragmentBases = getAlignedReadBases(read, 0);
            fragmentBases += getAlignedReadBases(mate, 0);

            int gapPosStart = readPosEnd < matePosStart ? readPosEnd : matePosEnd + 1;
            int gapPosEnd = readPosEnd < matePosStart ? matePosStart - 1 : readPosStart - 1;

            String gapRefBases = mConfig.RefGenome.getBaseString(read.getReferenceName(), gapPosStart, gapPosEnd);
            fragmentBases += gapRefBases;

            gcContent = GcCalcs.calcGcPercent(fragmentBases);
        }

        addFragmentData(fragPosStart, fragPosEnd, duplicateCount, gcContent, fragmentLength);
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

    private static double getReadGcContent(final SAMRecord read)
    {
        return GcCalcs.calcGcPercent(getAlignedReadBases(read, 0));
    }

    private void addFragmentData(int fragmentPosStart, int fragmentPosEnd, int duplicateCount, double rawGcPercent, int rawFragmentLength)
    {
        double gcPercent = roundGcPercent(rawGcPercent, mConfig.GcPercentUnits);
        int fragmentLength = roundFragmentLength(rawFragmentLength, mConfig.FragmentLengthUnits);

        int dupCountRounded = roundDuplicateCount(duplicateCount);

        mAllFragmentGcMap.add(fragmentLength, gcPercent, dupCountRounded);

        boolean matchesRegion = false;

        int fragBoundsStart = fragmentPosStart - TARGET_REGION_PROXIMITY;
        int fragBoundsEnd = fragmentPosEnd + TARGET_REGION_PROXIMITY;

        for(TargetRegionData targetRegion : mCurrentPartition.TargetRegions)
        {
            if(positionsOverlap(fragBoundsStart, fragBoundsEnd, targetRegion.start(), targetRegion.end()))
            {
                matchesRegion = true;

                if(mConfig.CaptureRegionCounts)
                    targetRegion.FragmentGcCounts.add(fragmentLength, gcPercent, dupCountRounded);
                else
                    break;
            }
        }

        if(matchesRegion)
            mTargetedFragmentGcMap.add(fragmentLength, gcPercent, dupCountRounded);
        else
            mNonTargetedFragmentGcMap.add(fragmentLength, gcPercent, dupCountRounded);
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

        int readIndex = 0;
        int position = read.getAlignmentStart();

        String alignedReadBases = "";

        for(CigarElement element : read.getCigar().getCigarElements())
        {
            switch(element.getOperator())
            {
                case S:
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
}
