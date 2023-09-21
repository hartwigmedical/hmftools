package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AnnotateTask implements Callable
{
    private final String mChromosome;
    private final AnnotateConfig mConfig;
    private final List<AnnotatedBedRecord> mBedRecords;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private final PerformanceCounter mPerfCounter;
    private int mRecordCounter;

    public AnnotateTask(final @NotNull String chromosome, @NotNull final AnnotateConfig config,
            @NotNull final List<AnnotatedBedRecord> bedRecords)
    {
        mChromosome = chromosome;
        mConfig = config;
        mBedRecords = bedRecords;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mPerfCounter = new PerformanceCounter(String.format("Chr %s", chromosome));
        mRecordCounter = 0;
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("Task for chromosome {} is starting.", mChromosome);

        RefGenomeCoordinates refGenomeCoords =
                mConfig.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        final int chromosomeLength = refGenomeCoords.length(stripChrPrefix(mChromosome));
        final ChrBaseRegion partition = new ChrBaseRegion(mChromosome, 1, chromosomeLength);

        mPerfCounter.start();
        mBamSlicer.slice(mSamReader, partition, this::processSamRecord);
        mPerfCounter.stop();

        mPerfCounter.logStats();
        MD_LOGGER.info("Task for chromosome is finished, {} reads processed", mChromosome, mRecordCounter);

        return (long) 0;
    }

    private void processSamRecord(@NotNull final SAMRecord read)
    {
        ++mRecordCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        for(var bedRecord : mBedRecords)
        {
            if(!bedRecord.getChromosome().equals(stripChrPrefix(mChromosome)))
            {
                MD_LOGGER.error("AnnotateTask for processing chromosome {} was passed a BED record associated to a different chromosome {}", stripChrPrefix(mChromosome), bedRecord.getChromosome());
                System.exit(1);
            }

            if(read.getAlignmentStart() >= bedRecord.getPosStart() && read.getAlignmentEnd() <= bedRecord.getPosEnd())
            {
                bedRecord.matchedRead(read);
            }
        }
    }
}
