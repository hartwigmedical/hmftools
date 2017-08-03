package com.hartwig.hmftools.purple.ratio;

import static com.hartwig.hmftools.common.cobalt.ReadCountFile.generateFilename;
import static com.hartwig.hmftools.common.cobalt.ReadCountFile.readFile;
import static com.hartwig.hmftools.common.copynumber.freec.FreecCPNFileLoader.normalReadCountLines;
import static com.hartwig.hmftools.common.copynumber.freec.FreecCPNFileLoader.tumorReadCountLines;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

class ReadCountSupplier {

    private static final Logger LOGGER = LogManager.getLogger(ReadCountSupplier.class);

    private final Multimap<String, ReadCount> tumorReadCount;
    private final Multimap<String, ReadCount> referenceReadCount;

    ReadCountSupplier(final CommonConfig config) throws IOException, HartwigException {

        final String cobaltTumorSampleFilename = generateFilename(config.cobaltDirectory(), config.tumorSample());
        final String cobaltReferenceSampleFilename = generateFilename(config.cobaltDirectory(), config.refSample());
        if (new File(cobaltTumorSampleFilename).exists()) {
            LOGGER.info("Loading cobalt read count data");
            tumorReadCount = readFile(cobaltTumorSampleFilename);
            referenceReadCount = readFile(cobaltReferenceSampleFilename);
        } else {
            LOGGER.info("Loading freec read count data");
            tumorReadCount = tumorReadCountLines(config.freecDirectory(), config.tumorSample());
            referenceReadCount = normalReadCountLines(config.freecDirectory(), config.refSample());
        }
    }

    Multimap<String, ReadCount> tumorReadCount() {
        return tumorReadCount;
    }

    Multimap<String, ReadCount> referenceReadCount() {
        return referenceReadCount;
    }

}
