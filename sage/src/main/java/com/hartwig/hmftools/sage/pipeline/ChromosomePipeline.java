package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.ReferenceData;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

import org.checkerframework.checker.units.qual.A;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ChromosomePipeline implements AutoCloseable
{
    private final String mChromosome;
    private final SageConfig mConfig;
    private final IndexedFastaSequenceFile mRefGenome;
    private final VcfWriter mVcfWriter;
    private final ChromosomePartition mPartition;
    private final Queue<RegionTask> mRegionTasks;
    private final RegionResults mRegionResults;

    private static final EnumSet<VariantTier> PANEL_ONLY_TIERS = EnumSet.of(VariantTier.HOTSPOT, VariantTier.PANEL);

    public ChromosomePipeline(
            final String chromosome, final SageConfig config,
            final ReferenceData refData, final Map<String,QualityRecalibrationMap> qualityRecalibrationMap,
            final Coverage coverage, final PhaseSetCounter phaseSetCounter, final VcfWriter vcfWriter)
    {
        mChromosome = chromosome;
        mConfig = config;
        mRefGenome = loadRefGenome(config.RefGenomeFile);
        mVcfWriter = vcfWriter;

        mRegionResults = new RegionResults(vcfWriter);

        final Chromosome chr = HumanChromosome.contains(chromosome)
                ? HumanChromosome.fromString(chromosome) : MitochondrialChromosome.fromString(chromosome);

        mPartition = new ChromosomePartition(config, mRefGenome);

        List<ChrBaseRegion> partitionedRegions = mPartition.partition(mChromosome);

        List<BaseRegion> chrPanel = refData.PanelWithHotspots.get(chr);
        List<VariantHotspot> chrHotspots = refData.Hotspots.get(chr);
        List<TranscriptData> chrTranscripts = refData.ChromosomeTranscripts.get(chromosome);
        List<BaseRegion> chrHighConfidence = refData.HighConfidence.get(chr);

        mRegionTasks= new ConcurrentLinkedQueue<>();
        // mRegionTasks = Lists.newArrayList();

        for(int i = 0; i < partitionedRegions.size(); ++i)
        {
            ChrBaseRegion region = partitionedRegions.get(i);

            List<BaseRegion> regionPanel = chrPanel != null ? chrPanel.stream()
                    .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                    : Lists.newArrayList();

            if(mConfig.PanelOnly && regionPanel.isEmpty())
                continue;

            List<VariantHotspot> regionHotspots = chrHotspots != null ? chrHotspots.stream()
                    .filter(x -> region.containsPosition(x.position())).collect(Collectors.toList()) : Lists.newArrayList();

            List<TranscriptData> regionsTranscripts = chrTranscripts != null ? chrTranscripts.stream()
                    .filter(x -> positionsOverlap(region.start(), region.end(), x.TransStart, x.TransEnd)).collect(Collectors.toList())
                    : Lists.newArrayList();

            List<BaseRegion> regionHighConfidence = chrHighConfidence != null ? chrHighConfidence.stream()
                    .filter(x -> positionsOverlap(region.start(), region.end(), x.start(), x.end())).collect(Collectors.toList())
                    : Lists.newArrayList();

            mRegionTasks.add(new RegionTask(
                    i, region, mRegionResults, config, mRefGenome, regionHotspots, regionPanel, regionsTranscripts,
                    regionHighConfidence, qualityRecalibrationMap, phaseSetCounter, coverage));
        }
    }

    public String chromosome()
    {
        return mChromosome;
    }

    public void process()
    {
        int regionCount = mRegionTasks.size();
        SG_LOGGER.info("chromosome({}) executing {} regions", mChromosome, regionCount);

        //final List<Callable> callableList = mRegionTasks.stream().collect(Collectors.toList());
        // TaskExecutor.executeTasks(callableList, mConfig.Threads);

        List<Thread> workers = new ArrayList<>();
        for (int i = 0; i < mConfig.Threads; ++i)
        {
            workers.add(new Thread(() -> {
                while(true)
                {
                    try
                    {
                        RegionTask task = mRegionTasks.remove();

                        if((mRegionTasks.size() % 100) == 0)
                        {
                            SG_LOGGER.debug("chromosome({}) regions completed({}) remaining({})",
                                    mChromosome, regionCount - mRegionTasks.size(), mRegionTasks.size());
                        }

                        task.call();
                    }
                    catch(NoSuchElementException e)
                    {
                        SG_LOGGER.trace("all tasks complete");
                        break;
                    }
                } }));
        }

        workers.forEach(x -> x.start());

        for(Thread worker : workers)
        {
            try
            {
                worker.join();
            }
            catch(InterruptedException e)
            {
                SG_LOGGER.error("task execution error: {}", e.toString());
                e.printStackTrace();
            }
        }

        SG_LOGGER.debug("chromosome({}) {} regions complete, processed {} reads, writing {} variants",
                mChromosome, regionCount, mRegionResults.totalReads(), mRegionResults.totalVariants());

        mVcfWriter.flushChromosome();

        mRegionResults.logPerfCounters();

        SG_LOGGER.info("chromosome({}) analysis complete", mChromosome);
    }

    @Override
    public void close() throws IOException
    {
        mRefGenome.close();
    }
}
