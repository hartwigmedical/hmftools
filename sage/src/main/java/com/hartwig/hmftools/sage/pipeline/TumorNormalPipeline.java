package com.hartwig.hmftools.sage.pipeline;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.NormalRefContextSupplier;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextCandidates;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.sam.SamSlicerFactoryChromImpl;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class TumorNormalPipeline {

    private static final Logger LOGGER = LogManager.getLogger(TumorNormalPipeline.class);

    private final String sample;
    private final String bamFile;
    private final String chromosome;
    private final SageConfig config;
    private final Executor executor;
    private final SamSlicerFactory samSlicerFactory;
    private final IndexedFastaSequenceFile refGenome;
    private final List<RefContext> result = Lists.newArrayList();
    private final RefContextCandidates candidates;

    private final List<CompletableFuture<List<RefContext>>> regions = Lists.newArrayList();

    public TumorNormalPipeline(@NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final List<GenomeRegion> panelRegions, final Executor executor,
            final RefContextCandidates candidates) {
        this.chromosome = chromosome;
        this.config = config;
        this.executor = executor;
        this.candidates = candidates;
        this.sample = config.reference();
        this.bamFile = config.referenceBam();
        this.samSlicerFactory = new SamSlicerFactoryChromImpl(config, panelRegions);
        this.refGenome = refGenome;

    }

    @NotNull
    public CompletableFuture<TumorNormalPipeline> submit(int minPosition, int maxPosition) {
        for (int i = 0; ; i++) {
            int start = 1 + i * config.regionSliceSize();
            int end = Math.min(start + config.regionSliceSize(), maxPosition);
            if (end >= minPosition) {
                start = Math.max(minPosition, start);
                addRegion(GenomeRegions.create(chromosome, start, end));
            }

            if (end >= maxPosition) {
                break;
            }
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(regions.toArray(new CompletableFuture[regions.size()]));
        return doneTumor.thenApply(aVoid -> {
            for (CompletableFuture<List<RefContext>> region : regions) {
                result.addAll(region.join());
            }
            return TumorNormalPipeline.this;
        });
    }

    @NotNull
    public List<RefContext> candidates() {
        return result;
    }

    private void addRegion(@NotNull final GenomeRegion region) {
        final RefSequence refSequence = new RefSequence(region, refGenome);
        final NormalRefContextSupplier supplier =
                new NormalRefContextSupplier(config, region, bamFile, refSequence, candidates, samSlicerFactory, refSequence);

        CompletableFuture<List<RefContext>> candidateFuture = CompletableFuture.supplyAsync(supplier, executor);
        regions.add(candidateFuture);
    }

}
