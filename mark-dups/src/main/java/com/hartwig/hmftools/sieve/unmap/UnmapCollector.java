package com.hartwig.hmftools.sieve.unmap;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.File;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class UnmapCollector implements Callable
{
    private final UnmapperConfig mConfig;
    private final ArrayBlockingQueue<ChrBaseRegion> mJobs;
    private final Set<PrimaryReadInfo> mReadsToDrop;

    private final RefGenomeSource mRefGenome;
    private final BamSlicer mBamSlicer;
    private final SamReader mSamReader;

    private ChrBaseRegion mCurrentRegion;

    private final PerformanceCounter mPerfCounter;
    private int mReadCounter;
    private int mPartitionCounter;

    public UnmapCollector(final UnmapperConfig config, final ArrayBlockingQueue<ChrBaseRegion> jobs, final Set<PrimaryReadInfo> readsToDrop)
    {
        mConfig = config;
        mJobs = jobs;
        mReadsToDrop = readsToDrop;

        mRefGenome = loadRefGenome(config.RefGenome);
        mBamSlicer = new BamSlicer(0, true, false, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mCurrentRegion = null;

        mPerfCounter = new PerformanceCounter("UnmapCollector task");
        mReadCounter = 0;
        mPartitionCounter = 0;
    }

    @Override
    public Long call()
    {
        while((mCurrentRegion = mJobs.poll()) != null)
        {
            ++mPartitionCounter;

            mPerfCounter.start();
            mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);
            mPerfCounter.stop();
        }

        mPerfCounter.logStats();
        UnmapperConfig.MD_LOGGER.info("UnmapCollector is finished, {} reads processed, {} regions processed", mReadCounter, mPartitionCounter);

        return (long) 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            UnmapperConfig.MD_LOGGER.info("UnmapCollector came across an unmapped read.");
            System.exit(1);
        }

        if(read.isSecondaryOrSupplementary())
        {
            UnmapperConfig.MD_LOGGER.info("UnmapCollector came across a secondary or supplementary read.");
            System.exit(1);
        }

        final ChrBaseRegion alignedRegion = new ChrBaseRegion(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
        if(!isHighFreqGRegion(alignedRegion))
        {
            return;
        }

        // TODO(m_cooper): Can I assume that all reads are paired?
        // TODO(m_cooper): Performance implications of various SamRecord methods.
        PrimaryReadInfo readInfo = new PrimaryReadInfo(read.getReadName(), read.getFirstOfPairFlag());
        Unmapper.registerReadForDropping(mReadsToDrop, readInfo);
    }

    // TODO(m_cooper): Performance implicates of random access.
    private boolean isHighFreqGRegion(final ChrBaseRegion region)
    {
        final String refBases = mRefGenome.getBaseString(region.Chromosome, region.start(), region.end());
        if(refBases.isEmpty())
        {
            MD_LOGGER.error("Requested empty base string from ref genome ({}) at {}:{}-{}.", mConfig.RefGenome, region.Chromosome, region.start(), region.end());
            System.exit(1);
        }

        int gCount = 0;
        for(int i = 0; i < refBases.length(); ++i)
        {
            if(refBases.charAt(i) == 'G')
            {
                ++gCount;
            }
        }

        // TODO(m_cooper): Make fraction configurable.
        return gCount >= 0.9 * refBases.length();
    }
}
