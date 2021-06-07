package com.hartwig.hmftools.sage.pipeline;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.function.Consumer;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneCoverage;
import com.hartwig.hmftools.sage.sam.SamSlicer;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CoveragePipeline
{

    private final SageConfig config;
    private final Coverage coverage;
    private final ChromosomePartition partition;
    private final ReferenceSequenceFile refGenome;
    private final ListMultimap<Chromosome, NamedBed> panel;
    private final ExecutorService executorService;

    public CoveragePipeline(final SageConfig config, final ReferenceSequenceFile refGenome, final ListMultimap<Chromosome, NamedBed> panel,
            final ExecutorService executorService)
    {
        this.executorService = executorService;

        final Set<String> samples = Sets.newHashSet();
        samples.addAll(config.reference());
        samples.addAll(config.tumor());

        this.config = config;
        this.partition = new ChromosomePartition(config, refGenome);
        this.panel = panel;
        this.refGenome = refGenome;
        this.coverage = new Coverage(samples, panel.values());
    }

    @NotNull
    public Coverage process() throws ExecutionException, InterruptedException
    {

        final List<Future<?>> futures = Lists.newArrayList();
        final SAMSequenceDictionary dictionary = refGenome.getSequenceDictionary();
        for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
        {
            final String contig = samSequenceRecord.getSequenceName();
            if(config.chromosomes().isEmpty() || config.chromosomes().contains(contig))
            {
                for(GenomeRegion region : partition.partition(contig))
                {

                    //                    for (int i = 0; i < config.tumor().size(); i++) {
                    //                        final String sample = config.tumor().get(i);
                    //                        final String sampleBam = config.tumorBam().get(i);
                    //                        futures.add(submit(sample, sampleBam, region));
                    //                    }

                    for(int i = 0; i < config.reference().size(); i++)
                    {
                        final String sample = config.reference().get(i);
                        final String sampleBam = config.referenceBam().get(i);
                        futures.add(submit(sample, sampleBam, region));
                    }
                }

            }
        }

        for(Future<?> future : futures)
        {
            future.get();
        }

        return coverage;
    }

    @NotNull
    private Future<?> submit(@NotNull final String sample, @NotNull final String bam, @NotNull final GenomeRegion bounds)
    {
        return executorService.submit(() ->
        {
            final List<GeneCoverage> geneCoverage = coverage.coverage(sample, bounds.chromosome());
            if(geneCoverage.isEmpty() || !HumanChromosome.contains(bounds.chromosome()))
            {
                return;
            }

            final Consumer<SAMRecord> consumer = record ->
            {
                final GenomeRegion alignment =
                        GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
                geneCoverage.forEach(x -> x.accept(alignment));
            };

            final HumanChromosome chrom = HumanChromosome.fromString(bounds.chromosome());
            final SamSlicer slicer = new SamSlicer(1, bounds, panel.get(chrom));
            try(final SamReader samReader = SamReaderFactory.makeDefault()
                    .validationStringency(config.validationStringency())
                    .referenceSource(new ReferenceSource(refGenome))
                    .open(new File(bam)))
            {
                slicer.slice(samReader, consumer);
            } catch(IOException e)
            {
                throw new CompletionException(e);
            }
        });
    }
}
