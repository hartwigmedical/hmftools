package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.google.common.util.concurrent.UncheckedExecutionException;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.BamSlicerFilter;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

// For performance reasons, this class is not decoupled from the threading
// We want to avoid locks, queues as much as possible
public class SampleBamProcessor
{
    private final BamSlicerFilter mBamSlicerFilter;
    private final SampleReadProcessor mSampleReadProcessor;
    private final List<ChrBaseRegion> mPartitions;

    public Collection<MicrosatelliteSiteAnalyser> getMicrosatelliteSiteAnalysers()
    {
        return mSampleReadProcessor.getMicrosatelliteSiteAnalysers();
    }

    public SampleBamProcessor(final JitterAnalyserConfig config, final BamSlicerFilter bamSlicerFilter,
            final SampleReadProcessor sampleReadProcessor)
    {
        mBamSlicerFilter = bamSlicerFilter;
        mSampleReadProcessor = sampleReadProcessor;

        // partition genome covered by MS analysers
        SortedSet<String> chromosomes = mSampleReadProcessor
                .getMicrosatelliteSiteAnalysers()
                .stream()
                .map(x -> x.refGenomeMicrosatellite.chromosome())
                .collect(Collectors.toCollection(() -> Sets.newTreeSet(Comparator.comparingInt(HumanChromosome::chromosomeRank)
                        .thenComparing(Function.identity()))));

        mPartitions = new ArrayList<>();
        for(String chromosome : chromosomes)
        {
            mPartitions.addAll(partitionChromosome(chromosome, config.RefGenVersion, config.PartitionSize));
        }
    }

    public void queryBam(final JitterAnalyserConfig config, ExecutorService executorService) throws InterruptedException
    {
        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }

        BamSlicer bamSlicer = new BamSlicer(mBamSlicerFilter);

        CompletableFuture<Void> bamSliceTasks = bamSlicer.queryAsync(new File(config.BamPath), readerFactory, mPartitions,
                false, executorService, mSampleReadProcessor::processRead);
        try
        {
            // wait for all to complete
            bamSliceTasks.get();
        }
        catch(ExecutionException e)
        {
            throw new UncheckedExecutionException(e);
        }
    }
}
