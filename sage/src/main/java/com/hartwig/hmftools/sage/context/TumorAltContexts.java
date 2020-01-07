package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.LocalPhaseSetAltContext;
import com.hartwig.hmftools.sage.phase.MnvCandidate;
import com.hartwig.hmftools.sage.phase.MnvCandidates;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.sam.SamSlicerFactoryChromImpl;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class TumorAltContexts {

    private static final Logger LOGGER = LogManager.getLogger(TumorAltContexts.class);

    private final String sample;
    private final String bamFile;
    private final String chromosome;
    private final SageConfig config;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final SamSlicerFactory samSlicerFactory;
    private final IndexedFastaSequenceFile refGenome;
    private final List<AltContext> mnv = Lists.newArrayList();
    private final List<AltContext> result = Lists.newArrayList();

    private final List<CompletableFuture<List<AltContext>>> regions = Lists.newArrayList();

    public TumorAltContexts(int sample, @NotNull final String chromosome, @NotNull final SageConfig config,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions, final Executor executor) {
        this.chromosome = chromosome;
        this.config = config;
        this.executor = executor;
        this.sample = config.tumor().get(sample);
        this.bamFile = config.tumorBam().get(sample);
        this.refGenome = refGenome;
        this.samSlicerFactory = new SamSlicerFactoryChromImpl(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
    }

    @NotNull
    public CompletableFuture<TumorAltContexts> submit() {
        return submit(1, refGenome.getSequence(chromosome).length());
    }

    @NotNull
    public CompletableFuture<TumorAltContexts> submit(int minPosition, int maxPosition) {

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

            final MnvCandidates mnvCandidates = new MnvCandidates(config, result::add, refGenome, panelRegions, hotspots);
            final LocalPhaseSetAltContext phase = new LocalPhaseSetAltContext(mnvCandidates);

            LOGGER.info("Mnv candidates {} chromosome {}", sample, chromosome);
            for (CompletableFuture<List<AltContext>> region : regions) {
                region.join().forEach(phase);
            }

            phase.flush();
            mnvCandidates.flush();

            MnvContextSupplier mnvContextSupplier =
                    new MnvContextSupplier(config, sample, bamFile, samSlicerFactory, hotspots, panelRegions, refGenome);
            for (MnvCandidate mnvCandidate : mnvCandidates.mvnCandidates()) {
                mnv.addAll(mnvContextSupplier.get(mnvCandidate));
            }

            return TumorAltContexts.this;
        });
    }

    private void addRegion(@NotNull final GenomeRegion region) {
        final RefSequence refSequence = new RefSequence(region, refGenome);
        final AltContextSupplier supplier =
                new AltContextSupplier(config, sample, region, bamFile, refSequence, samSlicerFactory, refSequence, hotspots, panelRegions);

        CompletableFuture<List<AltContext>> candidateFuture = CompletableFuture.supplyAsync(supplier, executor);
        regions.add(candidateFuture);
    }

}
