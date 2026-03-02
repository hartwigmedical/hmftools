package com.hartwig.hmftools.bamtools.utils;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.utils.BamUmiAnalyser.writeReadInfo;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionUmiAnalyser extends Thread
{
    private final UmiAnalyserConfig mConfig;
    private final TaskQueue mRegions;
    private final Set<String> mKnownUmiSet;
    private final BufferedWriter mWriter;

    private final SamReader mUmiSamReader;
    private final SamReader mRawSamReader;

    private ChrBaseRegion mCurrentRegion;

    private final Map<String,SAMRecord> mReadCache;

    private final int mMaxKnownUmiLength;
    private long mUmiReadCount;
    private long mRawReadCount;

    private final boolean mLogReadIds;

    public PartitionUmiAnalyser(
            final UmiAnalyserConfig config, final TaskQueue regions, final Set<String> knownUmiSet, final BufferedWriter writer)
    {
        mConfig = config;
        mRegions = regions;
        mKnownUmiSet = knownUmiSet;
        mMaxKnownUmiLength = mKnownUmiSet.stream().mapToInt(x -> x.length()).max().orElse(0);
        mWriter = writer;

        mCurrentRegion = null;

        mUmiSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.UmiBamFile));
        mRawSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.RawBamFile));

        mReadCache = Maps.newHashMap();
        mUmiReadCount = 0;
        mRawReadCount = 0;

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion region = (ChrBaseRegion)mRegions.removeItem();

                processRegion(region);
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void processRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;

        BamSlicer bamSlicer = new BamSlicer(0, true, false, false);
        bamSlicer.slice(mUmiSamReader, mCurrentRegion, this::processUmiRecord);

        bamSlicer.slice(mRawSamReader, mCurrentRegion, this::processRawRecord);

        BT_LOGGER.debug("region({}) complete, reads(umi={}, raw={})", mCurrentRegion, mUmiReadCount, mRawReadCount);

        if(!mReadCache.isEmpty())
        {
            BT_LOGGER.warn("region({}) purged {} unmatched reads", mCurrentRegion, mReadCache.size());
        }
    }

    private void processUmiRecord(final SAMRecord umiRead)
    {
        if(!mCurrentRegion.containsPosition(umiRead.getAlignmentStart()))
            return;

        ++mUmiReadCount;

        if(mLogReadIds && mConfig.LogReadIds.contains(umiRead.getReadName()))
        {
            BT_LOGGER.debug("specific read: {}", readToString(umiRead));
        }

        String readId = umiRead.getReadName();
        int delimIndex = readId.lastIndexOf(UMI_DELIM);
        String noUmiReadId = readId.substring(0, delimIndex);
        boolean firstInPair = umiRead.getFirstOfPairFlag();
        String cachedReadId = format("%s_%d", noUmiReadId, firstInPair ? 1 : 2);

        mReadCache.put(cachedReadId, umiRead);
    }

    private static final String UMI_DELIM = ":";

    private void processRawRecord(final SAMRecord rawRead)
    {
        if(!mCurrentRegion.containsPosition(rawRead.getAlignmentStart()))
            return;

        ++mRawReadCount;

        if(mLogReadIds && mConfig.LogReadIds.contains(rawRead.getReadName()))
        {
            BT_LOGGER.debug("specific read: {}", readToString(rawRead));
        }

        // find matching read
        String readId = rawRead.getReadName();
        boolean firstInPair = rawRead.getFirstOfPairFlag();
        String cachedReadId = format("%s_%d", readId, firstInPair ? 1 : 2);

        SAMRecord umiRead = mReadCache.remove(cachedReadId);

        if(umiRead == null)
        {
            BT_LOGGER.debug("unmatched read: {}", readToString(rawRead));
            return;
        }

        // find the extracted bases
        String umiReadBases = umiRead.getReadString();
        String rawReadBases = rawRead.getReadString();

        int baseMatchIndex = rawReadBases.indexOf(umiReadBases);

        if(baseMatchIndex < 0)
        {
            BT_LOGGER.debug("read({}) unmatched bases, details({})", readId, readToString(rawRead));
            BT_LOGGER.debug("umi read({}) bases: {}", umiRead.getReadName(), umiReadBases);
            BT_LOGGER.debug("raw read({}) bases: {}", rawRead.getReadName(), rawReadBases);
            return;
        }

        // the UMI is stripped from the 5' end, so check this
        boolean expectUmiAtStart = !rawRead.getReadNegativeStrandFlag();

        int readBaseLengthDiff = rawReadBases.length() - umiReadBases.length();

        int leftStrippedBases = baseMatchIndex;
        int rightStrippedBases = max(readBaseLengthDiff - baseMatchIndex, 0);

        if((expectUmiAtStart && leftStrippedBases == 0) || (!expectUmiAtStart && rightStrippedBases == 0))
        {
            BT_LOGGER.warn("read({}) matchIndex({}) at wrong end details({})", readId, baseMatchIndex, readToString(rawRead));
        }

        int umiLength = 0;
        String umiSequence = "";
        int otherStrippedLength;

        if(expectUmiAtStart)
        {
            umiLength = leftStrippedBases;
            umiSequence = rawReadBases.substring(0, leftStrippedBases);
            otherStrippedLength = rightStrippedBases;
        }
        else
        {
            umiLength = rightStrippedBases;
            umiSequence = rawReadBases.substring(rawReadBases.length() - umiLength);
            otherStrippedLength = leftStrippedBases;
        }

        if(mMaxKnownUmiLength > 0 && umiLength > mMaxKnownUmiLength)
        {

        }

        int delimIndex = umiRead.getReadName().lastIndexOf(UMI_DELIM);
        String umiReadIdTag = umiRead.getReadName().substring(delimIndex + 1);

        boolean sequenceMatched = hasKnownUmiMatch(umiSequence);

        String reverseUmiSequence = Nucleotides.reverseComplementBases(umiSequence);
        boolean revRequenceMatched = hasKnownUmiMatch(reverseUmiSequence);

        writeReadInfo(
                mWriter, rawRead, umiReadIdTag, umiSequence, sequenceMatched, reverseUmiSequence, revRequenceMatched, otherStrippedLength);
    }

    private boolean hasKnownUmiMatch(final String umi)
    {
        if(mKnownUmiSet.contains(umi))
            return true;

        for(String knownUmi : mKnownUmiSet)
        {
            if(!exceedsUmiDiff(umi, knownUmi, 1))
                return true;
        }

        return false;
    }

    private static boolean exceedsUmiDiff(final String first, final String second, int permittedDiff)
    {
        if(first.length() != second.length())
            return true;

        short diffs = 0;
        for(short i = 0; i < first.length(); ++i)
        {
            if(first.charAt(i) != second.charAt(i))
            {
                ++diffs;

                if(diffs > permittedDiff)
                    return true;
            }
        }

        return false;
    }
}
