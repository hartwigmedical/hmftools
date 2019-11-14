package com.hartwig.hmftools.patientreporter.variants.somatic;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;

import org.jetbrains.annotations.NotNull;

final class MicrosatelliteAnalyzer {

    private MicrosatelliteAnalyzer() {
    }

    static double determineMicrosatelliteIndelsPerMb(@NotNull List<SomaticVariant> variants) {
        final MicrosatelliteIndels indelCount = new MicrosatelliteIndels();
        variants.forEach(indelCount::accept);
        return indelCount.microsatelliteIndelsPerMb();
    }
}
