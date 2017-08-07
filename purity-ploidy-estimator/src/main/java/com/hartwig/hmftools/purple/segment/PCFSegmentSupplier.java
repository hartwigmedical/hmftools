package com.hartwig.hmftools.purple.segment;

import java.io.File;
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
import com.hartwig.hmftools.common.purple.baf.TumorBAFFile;
import com.hartwig.hmftools.common.purple.pcf.ImmutablePCFRegion;
import com.hartwig.hmftools.common.purple.pcf.PCFFile;
import com.hartwig.hmftools.common.purple.pcf.PCFRegion;
import com.hartwig.hmftools.common.purple.ratio.ReadRatioFile;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.purple.PurityPloidyEstimateApplication;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.rosuda.REngine.REXPMismatchException;
import org.rosuda.REngine.Rserve.RConnection;
import org.rosuda.REngine.Rserve.RserveException;

public class PCFSegmentSupplier implements Supplier<List<GenomeRegion>> {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    private final List<GenomeRegion> segments = Lists.newArrayList();

    public PCFSegmentSupplier(@NotNull final ExecutorService executorService, @NotNull final CommonConfig config,
            Map<String, ChromosomeLength> chromosomeLengths)
            throws RserveException, REXPMismatchException, IOException, ExecutionException, InterruptedException {

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
            @NotNull final String sample) throws RserveException, IOException {
        final String ratioFile = ReadRatioFile.generateFilename(config.outputDirectory(), sample);
        final String pcfFile = PCFFile.generateRatioFilename(config.outputDirectory(), sample);
        if (config.forceRecalculate() || !new File(pcfFile).exists()) {
            LOGGER.info("Connecting to R server to generate {} segmentation", description);

            RConnection c = new RConnection();

            c.eval("library(copynumber)");
            c.eval("ratio <- read.table(\"" + ratioFile + "\", header=TRUE)");
            c.eval("ratio <- ratio[ratio$Ratio>=0,]");
            c.eval("ratio$Ratio[ratio$Ratio<0.001] <- 0.001 ");
            c.eval("ratio$S1 = log2(ratio$Ratio)");
            c.eval("ratio <- ratio[!is.nan(ratio$S1),]");
            c.eval("ratio <- ratio[,c(\"Chromosome\",\"Position\",\"S1\")]");
            c.eval("ratio.seg<-pcf(ratio,verbose=FALSE,gamma=100, kmin=1,save.res = TRUE, file.names = c(\"" + pcfFile + "1\", \"" + pcfFile
                    + "\"))");
        }

        LOGGER.info("Loading {} PCF segments from {}", description, pcfFile);
        return PCFFile.read(config.windowSize(), pcfFile);
    }

    @NotNull
    private Multimap<String, PCFRegion> bafSegmentation(final @NotNull CommonConfig config, @NotNull final String sample)
            throws RserveException, IOException {
        final String bafFile = TumorBAFFile.generateFilename(config.outputDirectory(), sample);
        final String pcfFile = PCFFile.generateBAFFilename(config.outputDirectory(), sample);
        if (config.forceRecalculate() || !new File(pcfFile).exists()) {
            LOGGER.info("Connecting to R server to generate baf segmentation");

            RConnection c = new RConnection();

            c.eval("library(copynumber)");
            c.eval("baf <- read.table(\"" + bafFile + "\", header=TRUE)");
            c.eval("baf <- baf[,c(\"Chromosome\",\"Position\",\"mBAF\")]");
            c.eval("baf.seg<-pcf(baf,verbose=FALSE,gamma=100, kmin=1,save.res = TRUE, file.names = c(\"" + pcfFile + "1\", \"" + pcfFile
                    + "\"))");
        }

        LOGGER.info("Loading baf PCF segments from {}", pcfFile);
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
