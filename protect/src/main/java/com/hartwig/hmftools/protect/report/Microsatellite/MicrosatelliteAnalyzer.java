package com.hartwig.hmftools.protect.report.Microsatellite;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;

import org.jetbrains.annotations.NotNull;

public final class MicrosatelliteAnalyzer {

    private MicrosatelliteAnalyzer() {
    }

    public static double determineMicrosatelliteIndelsPerMb(@NotNull List<SomaticVariant> variants) {
        final MicrosatelliteIndels indelCount = new MicrosatelliteIndels();
        variants.forEach(indelCount::accept);
        return indelCount.microsatelliteIndelsPerMb();
    }
}
