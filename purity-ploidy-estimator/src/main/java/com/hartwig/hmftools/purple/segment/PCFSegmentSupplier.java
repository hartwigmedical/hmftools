package com.hartwig.hmftools.purple.segment;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.ChromosomeLength;
import com.hartwig.hmftools.common.pcf.ImmutablePCFRegion;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.pcf.PCFRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.purple.PurityPloidyEstimateApplication;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PCFSegmentSupplier implements Supplier<List<GenomeRegion>> {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    private final List<GenomeRegion> segments = Lists.newArrayList();

    public PCFSegmentSupplier(@NotNull final ExecutorService executorService, @NotNull final CommonConfig config,
            Map<String, ChromosomeLength> chromosomeLengths) throws IOException, ExecutionException, InterruptedException {

        final Future<Multimap<String, PCFRegion>> tumorRegionsFuture =
                executorService.submit(callable(config, "tumor ratio", config.tumorSample()));
        final Future<Multimap<String, PCFRegion>> referenceRegionsFuture =
                executorService.submit(callable(config, "reference ratio", config.refSample()));

        Multimap<String, PCFRegion> tumorRegions = tumorRegionsFuture.get();
        Multimap<String, PCFRegion> referenceRegions = referenceRegionsFuture.get();

        for (String chromosome : referenceRegions.keySet()) {
            long chromosomeEnd = chromosomeLengths.get(chromosome).position();
            final List<GenomeRegion> tumor = Lists.newArrayList(tumorRegions.get(chromosome));
            final List<GenomeRegion> reference = extend(chromosome, chromosomeEnd, referenceRegions.get(chromosome));

            segments.addAll(SegmentMerge.merge(reference, tumor));
        }

        Collections.sort(segments);
    }

    private Callable<Multimap<String, PCFRegion>> callable(final @NotNull CommonConfig config, @NotNull final String description,
            @NotNull final String sample) {
        return () -> ratioSegmentation(config, description, sample);
    }

    @NotNull
    private Multimap<String, PCFRegion> ratioSegmentation(final @NotNull CommonConfig config, @NotNull final String description,
            @NotNull final String sample) throws IOException, InterruptedException {
        final String pcfFile = PCFFile.generateRatioFilename(config.cobaltDirectory(), sample);

        LOGGER.info("Loading {} PCF segments from {}", description, pcfFile);
        return PCFFile.read(config.windowSize(), pcfFile);
    }

    @Override
    public List<GenomeRegion> get() {
        return segments;
    }

    @NotNull
    private List<GenomeRegion> extend(final String chromosome, final long chromosomeEnd, Collection<PCFRegion> regions) {
        final List<GenomeRegion> result = Lists.newArrayList();
        long start = 1;
        long end = 1;
        for (GenomeRegion genomeRegion : regions) {
            if (genomeRegion.start() > start) {
                result.add(create(chromosome, start, genomeRegion.start() - 1));
            }

            end = genomeRegion.end();
            start = end + 1;
            result.add(create(chromosome, genomeRegion.start(), end));
        }

        if (end < chromosomeEnd) {
            result.add(create(chromosome, start, chromosomeEnd));
        }

        return result;
    }

    private GenomeRegion create(String chromosome, long start, long end) {
        return ImmutablePCFRegion.builder().chromosome(chromosome).start(start).end(end).build();
    }
}
