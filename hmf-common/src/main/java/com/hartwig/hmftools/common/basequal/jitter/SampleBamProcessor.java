package com.hartwig.hmftools.common.basequal.jitter;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.util.concurrent.UncheckedExecutionException;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;

// For performance reasons, this class is not decoupled from the threading
// We want to avoid locks, queues as much as possible
public class SampleBamProcessor
{
    private final List<RefGenomeMicrosatellite> mRefGenomeMicrosatellites;

    // for speed reasons we need to consolidate the chr base regions into bigger chunks
    private final Multimap<ChrBaseRegion, MicrosatelliteSiteAnalyser> mMicrosatelliteSiteAnalysers = ArrayListMultimap.create();

    public Collection<MicrosatelliteSiteAnalyser> getMicrosatelliteSiteAnalysers() { return mMicrosatelliteSiteAnalysers.values(); }

    public SampleBamProcessor(List<RefGenomeMicrosatellite> refGenomeMicrosatellites)
    {
        mRefGenomeMicrosatellites = refGenomeMicrosatellites;
        partitionGenome();
    }

    private void partitionGenome()
    {
        final int PARTITION_SIZE = 1_000;

        mMicrosatelliteSiteAnalysers.clear();

        ImmutableListMultimap<String, RefGenomeMicrosatellite> chromosomeMsSites = Multimaps.index(mRefGenomeMicrosatellites, RefGenomeMicrosatellite::chromosome);

        List<MicrosatelliteSiteAnalyser> regionAnalysers = new ArrayList<>();

        for(String chromosome : chromosomeMsSites.keySet())
        {
            ChrBaseRegion currentRegion = null;

            List<RefGenomeMicrosatellite> refGenomeMicrosatellites = new ArrayList<>(chromosomeMsSites.get(chromosome));
            refGenomeMicrosatellites.sort(Comparator.comparing(RefGenomeMicrosatellite::referenceStart));

            for(RefGenomeMicrosatellite refGenomeMicrosatellite : refGenomeMicrosatellites)
            {
                if(currentRegion == null || currentRegion.baseLength() >= PARTITION_SIZE)
                {
                    if(currentRegion != null)
                    {
                        mMicrosatelliteSiteAnalysers.putAll(currentRegion, regionAnalysers);
                    }

                    currentRegion = refGenomeMicrosatellite.genomeRegion.clone();
                    regionAnalysers.clear();
                }
                currentRegion.setEnd(refGenomeMicrosatellite.genomeRegion.end());
                regionAnalysers.add(new MicrosatelliteSiteAnalyser(refGenomeMicrosatellite));
            }

            // final one
            if(currentRegion != null)
            {
                mMicrosatelliteSiteAnalysers.putAll(currentRegion, regionAnalysers);
            }
        }
    }

    public void queryBam(final JitterAnalyserConfig config, ExecutorService executorService) throws InterruptedException
    {
        SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(config.BamStringency);
        if(config.RefGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSource(new ReferenceSource(new File(config.RefGenomeFile)));
        }

        Collection<ChrBaseRegion> partitions = mMicrosatelliteSiteAnalysers.keySet().stream().sorted().collect(Collectors.toList());

        BamSlicer bamSlicer = new BamSlicer(config.MinMappingQuality, false, false, false);
        CompletableFuture<Void> bamSliceTasks = bamSlicer.queryAsync(new File(config.BamPath), readerFactory, partitions,
                false, executorService, this::processRead, this::regionComplete);
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

    public void processRead(SAMRecord read, ChrBaseRegion baseRegion)
    {
        if(read.getReadUnmappedFlag())
        {
            return;
        }

        Collection<MicrosatelliteSiteAnalyser> microsatelliteSiteAnalysers = this.mMicrosatelliteSiteAnalysers.get(baseRegion);

        Validate.isTrue(!microsatelliteSiteAnalysers.isEmpty());

        int readAlignmentStart = read.getAlignmentStart();
        int readAlignmentEnd = read.getAlignmentEnd();

        for(MicrosatelliteSiteAnalyser analyser : microsatelliteSiteAnalysers)
        {
            if(BaseRegion.positionsWithin(analyser.refGenomeMicrosatellite.referenceStart(), analyser.refGenomeMicrosatellite.referenceEnd(),
                    readAlignmentStart, readAlignmentEnd))
            {
                analyser.addReadToStats(read);
            }
        }
    }

    public void regionComplete(ChrBaseRegion baseRegion)
    {
    }
}
