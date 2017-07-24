package com.hartwig.hmftools.purple.ratio;

import static com.hartwig.hmftools.common.copynumber.freec.FreecCPNFileLoader.normalReadCountLines;
import static com.hartwig.hmftools.common.copynumber.freec.FreecCPNFileLoader.tumorReadCountLines;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.ratio.GCContent;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatios;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatiosBuilder;
import com.hartwig.hmftools.common.purple.ratio.ReadCount;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;
import com.hartwig.hmftools.common.purple.ratio.ReadRatioFile;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ReadCountRatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(ReadCountRatioSupplier.class);

    private final Multimap<String, ReadRatio> referenceRatios;
    private final Multimap<String, ReadRatio> tumorRatios;

    public ReadCountRatioSupplier(final CommonConfig config, final Multimap<String, GCContent> gcContent)
            throws IOException, HartwigException {

        LOGGER.info("Loading read count data");
        final Multimap<String, ReadCount> tumorReadCount = tumorReadCountLines(config.freecDirectory(), config.tumorSample());
        final Multimap<String, ReadCount> normalReadCount = normalReadCountLines(config.freecDirectory(), config.refSample());

        LOGGER.info("Generating gc normalized read ratios");
        final NormalizedRatiosBuilder normalRatiosBuilder = new NormalizedRatiosBuilder();
        final NormalizedRatiosBuilder tumorRatiosBuilder = new NormalizedRatiosBuilder();
        for (String chromosome : normalReadCount.keySet()) {
            List<GCContent> chromosomeGCContent = (List<GCContent>) gcContent.get(chromosome);
            List<ReadCount> chromosomeNormalReadCount = (List<ReadCount>) normalReadCount.get(chromosome);
            List<ReadCount> chromosomeTumorReadCount = (List<ReadCount>) tumorReadCount.get(chromosome);
            for (int i = 0; i < chromosomeNormalReadCount.size(); i++) {
                GCContent windowGCContent = chromosomeGCContent.get(i);
                ReadCount windowNormalCount = chromosomeNormalReadCount.get(i);
                ReadCount windowTumorCount = chromosomeTumorReadCount.get(i);

                normalRatiosBuilder.addPosition(windowGCContent, windowNormalCount);
                tumorRatiosBuilder.addPosition(windowGCContent, windowTumorCount);
            }
        }

        NormalizedRatios normalizedReferenceRatios = normalRatiosBuilder.build();
        NormalizedRatios normalizedTumorRatios = tumorRatiosBuilder.build();

        referenceRatios = normalizedReferenceRatios.normalisedRatios();
        tumorRatios = normalizedTumorRatios.normalisedRatios();

        LOGGER.info("Persisting gc normalized read ratios to file");
        ReadRatioFile.write(config.outputDirectory(), config.refSample(), referenceRatios);
        ReadRatioFile.write(config.outputDirectory(), config.tumorSample(), tumorRatios);
    }

    public Multimap<String, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    public Multimap<String, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

}
