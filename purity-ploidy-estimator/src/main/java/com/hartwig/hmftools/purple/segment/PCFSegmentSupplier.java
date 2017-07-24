package com.hartwig.hmftools.purple.segment;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.pcf.ImmutablePCFRegion;
import com.hartwig.hmftools.common.purple.pcf.PCFFile;
import com.hartwig.hmftools.common.purple.pcf.PCFRegion;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;
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
    private static final int WINDOW_SIZE = 1000;

    private final List<GenomeRegion> segments = Lists.newArrayList();

    public PCFSegmentSupplier(@NotNull final CommonConfig config, Multimap<String, ReadRatio> ratios)
            throws RserveException, REXPMismatchException, IOException {

        Multimap<String, PCFRegion> referenceRegions = getStringPCFRegionMultimap(config, config.refSample());
        Multimap<String, PCFRegion> tumorRegions = getStringPCFRegionMultimap(config, config.tumorSample());

        for (String chromosome : referenceRegions.keySet()) {
            long chromosomeEnd = ratios.get(chromosome).size() * WINDOW_SIZE;
            segments.addAll(extend(chromosome, chromosomeEnd, referenceRegions.get(chromosome)));
        }

        Collections.sort(segments);
    }

    @NotNull
    private Multimap<String, PCFRegion> getStringPCFRegionMultimap(final @NotNull CommonConfig config, String sample)
            throws RserveException, IOException {
        final String ratioFile = ReadRatioFile.generateFilename(config.outputDirectory(), sample);
        final String pcfFile = PCFFile.generateFilename(config.outputDirectory(), sample);
        if (!new File(pcfFile).exists()) {
            LOGGER.info("Connecting to R server to compute segmentation");

            RConnection c = new RConnection();
            String cmdLoadData = "ratio <- read.table(\"" + ratioFile + "\", header=TRUE)";
            String cmdCleanData =
                    "ratio <- ratio[ratio$Ratio>0,] \n ratio$S1 = log2(ratio$Ratio) \n ratio <- ratio[!is.nan(ratio$S1),] \n ratio <- ratio[,c(\"Chromosome\",\"Position\",\"S1\")]";
            String cmdProcessData =
                    "ratio.seg<-pcf(ratio,verbose=FALSE,gamma=100, kmin=2,save.res = TRUE, file.names = c(\"" + pcfFile + "1\", \""
                            + pcfFile + "\"))";

            c.eval("library(copynumber)");
            c.eval(cmdLoadData);
            c.eval(cmdCleanData);
            c.eval(cmdProcessData);
        }

        LOGGER.info("Loading PCF segments from {}", pcfFile);
        return PCFFile.read(pcfFile);
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

            end = genomeRegion.end() + WINDOW_SIZE - 1;
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
