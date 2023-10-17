package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
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
    private final ArrayBlockingQueue<ChrBaseRegion> mJobs;
    private final MaskedRegions mMaskedRegions;
    private final EnsemblDataCache mEnsemblDataCache;
    private final BufferedWriter mOutputWriter;

    private final RefGenomeSource mRefGenome;
    private final RefGenomeCoordinates mRefGenomeCoords;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mConsumerPerfCounter;
    private final PerformanceCounter mJobPerfCounter;

    private ChrBaseRegion mCurrentRegion;
    private HighDepthCounts mCurrentCounts;
    private long mReadCounter;
    private long mRegionCounter;

    public HighDepthCountConsumer(
            final AnnotateConfig config,
            final ArrayBlockingQueue<ChrBaseRegion> jobs,
            final MaskedRegions maskedRegions,
            final EnsemblDataCache ensemblDataCache,
            final BufferedWriter outputWriter)
    {
        mConfig = config;
        mJobs = jobs;
        mMaskedRegions = maskedRegions;
        mEnsemblDataCache = ensemblDataCache;
        mOutputWriter = outputWriter;

        mRefGenome = loadRefGenome(config.RefGenome);
        mRefGenomeCoords = mConfig.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        mBamSlicer = new BamSlicer(0, config.KeepDuplicates, true, false);
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

            mJobPerfCounter.start();
            RefGenomeRegionAnnotations refAnnotations = annotateRegionFromRefGenome();
            mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);
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
        String refBases = mRefGenome.getBaseString(
                mCurrentRegion.Chromosome,
                mCurrentRegion.start(),
                mCurrentRegion.end());

        if(refBases.isEmpty())
        {
            MD_LOGGER.error("Requested empty base string from ref genome ({}) at {}:{}-{}.", mConfig.RefGenome, mCurrentRegion.Chromosome, mCurrentRegion.start(), mCurrentRegion.end());
            System.exit(1);
        }

        // Base counts.
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
                        mCurrentRegion.Chromosome,
                        mCurrentRegion.start(),
                        mCurrentRegion.end());
                System.exit(1);
            }
        }

        // Distance to centromere.
        final int centromere = mRefGenomeCoords.centromere(mCurrentRegion.Chromosome);
        int distToCentromere = 0;
        if(mCurrentRegion.end() <= centromere)
        {
            distToCentromere = centromere - mCurrentRegion.end();
        }
        else if(centromere <= mCurrentRegion.start())
        {
            distToCentromere = mCurrentRegion.start() - centromere;
        }

        // Distance to telomere.
        final int chrLength = mRefGenomeCoords.length(mCurrentRegion.Chromosome);
        final int distToTelomere = Math.min(mCurrentRegion.start() - 1, chrLength - mCurrentRegion.end());

        // Find gene overlaps.
        final List<GeneData> genes =
                mEnsemblDataCache.findGenesByExonRegionOverlap(mCurrentRegion.Chromosome, mCurrentRegion.start(), mCurrentRegion.end());
        final String genesStr = genes.stream().map(geneData -> geneData.GeneName).collect(Collectors.joining(","));

        return new RefGenomeRegionAnnotations(
                distToCentromere,
                distToTelomere,
                mMaskedRegions.distance(mCurrentRegion),
                genesStr,
                1.0 * aCount / refBases.length(),
                1.0 * tCount / refBases.length(),
                1.0 * gCount / refBases.length(),
                1.0 * cCount / refBases.length());
    }
}
