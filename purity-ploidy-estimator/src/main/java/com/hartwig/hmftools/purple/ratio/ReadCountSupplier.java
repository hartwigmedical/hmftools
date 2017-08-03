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

class ReadCountSupplier {

    private final Multimap<String, ReadCount> tumorReadCount;
    private final Multimap<String, ReadCount> referenceReadCount;

    ReadCountSupplier(final CommonConfig config) throws IOException, HartwigException {

        final String cobaltTumorSampleFilename = generateFilename(config.cobaltDirectory(), config.tumorSample());
        final String cobaltReferenceSampleFilename = generateFilename(config.cobaltDirectory(), config.refSample());
        if (new File(cobaltTumorSampleFilename).exists()) {
            tumorReadCount = readFile(cobaltTumorSampleFilename);
            referenceReadCount = readFile(cobaltReferenceSampleFilename);
        } else {
            tumorReadCount = tumorReadCountLines(config.freecDirectory(), config.tumorSample());
            referenceReadCount = normalReadCountLines(config.freecDirectory(), config.refSample());
        }
    }

    public Multimap<String, ReadCount> tumorReadCount() {
        return tumorReadCount;
    }

    public Multimap<String, ReadCount> referenceReadCount() {
        return referenceReadCount;
    }

}
