package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class HighDepthCountConsumer implements Callable
{
    private final AnnotateConfig mConfig;
    private final ArrayBlockingQueue<HighDepthRegion> mJobs;
    private final BufferedWriter mOutputWriter;

    private final RefGenomeSource mRefGenome;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mConsumerPerfCounter;
    private final PerformanceCounter mJobPerfCounter;

    private HighDepthRegion mCurrentRegion;
    private HighDepthCounts mCurrentCounts;
    private long mReadCounter;
    private long mRegionCounter;

    public HighDepthCountConsumer(final AnnotateConfig config, final ArrayBlockingQueue<HighDepthRegion> jobs,
            final BufferedWriter outputWriter)
    {
        mConfig = config;
        mJobs = jobs;
        mOutputWriter = outputWriter;

        mRefGenome = loadRefGenome(config.RefGenome);
        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mJobPerfCounter = new PerformanceCounter("HighDepthCountConsumer Jobs");
        mConsumerPerfCounter = new PerformanceCounter("HighDepthCountConsumer Total");
        mReadCounter = 0;
        mRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        mConsumerPerfCounter.start();

        while((mCurrentRegion = mJobs.poll()) != null)
        {
            ++mRegionCounter;
            mCurrentCounts = new HighDepthCounts();
            final ChrBaseRegion chrBaseRegion =
                    new ChrBaseRegion(mCurrentRegion.getChromosome(), mCurrentRegion.getPosStart(), mCurrentRegion.getPosEnd());

            mJobPerfCounter.start();
            RefGenomeRegionAnnotations refAnnotations = annotateRegionFromRefGenome();
            mBamSlicer.slice(mSamReader, chrBaseRegion, this::processSamRecord);
            mJobPerfCounter.stop();

            Annotate.writeRecord(mOutputWriter, mCurrentRegion, refAnnotations, mCurrentCounts);
        }

        mConsumerPerfCounter.stop();

        mJobPerfCounter.logStats();
        mConsumerPerfCounter.logStats();
        MD_LOGGER.info("HighDepthCountConsumer is finished, {} reads processed, {} regions processed", mReadCounter, mRegionCounter);

        return (long) 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        mCurrentCounts.matchedRead(read);
    }

    // TODO(m_cooper): This is probably implemented somewhere.
    private RefGenomeRegionAnnotations annotateRegionFromRefGenome()
    {
        String refBases =
                mRefGenome.getBaseString(mCurrentRegion.getChromosome(), mCurrentRegion.getPosStart(), mCurrentRegion.getPosEnd());
        if(refBases.isEmpty())
        {
            MD_LOGGER.error("Requested empty base string from ref genome ({}) at {}:{}-{}.", mConfig.RefGenome, mCurrentRegion.getChromosome(), mCurrentRegion.getPosStart(), mCurrentRegion.getPosEnd());
            System.exit(1);
        }

        int aCount = 0;
        int tCount = 0;
        int gCount = 0;
        int cCount = 0;
        for(int i = 0; i < refBases.length(); i++)
        {
            if(refBases.charAt(i) == 'A')
            {
                aCount++;
            }
            else if(refBases.charAt(i) == 'T')
            {
                tCount++;
            }
            else if(refBases.charAt(i) == 'G')
            {
                gCount++;
            }
            else if(refBases.charAt(i) == 'C')
            {
                cCount++;
            }
            else if(refBases.charAt(i) == 'N')
            {
            }
            else
            {
                MD_LOGGER.error(
                        "Found an unknown base ({}) in the ref genome ({}) at {}:{}-{}.",
                        refBases.charAt(i),
                        mConfig.RefGenome,
                        mCurrentRegion.getChromosome(),
                        mCurrentRegion.getPosStart(),
                        mCurrentRegion.getPosEnd());
                System.exit(1);
            }
        }

        return new RefGenomeRegionAnnotations(
                1.0 * aCount / refBases.length(),
                1.0 * tCount / refBases.length(),
                1.0 * gCount / refBases.length(),
                1.0 * cCount / refBases.length());
    }
}
