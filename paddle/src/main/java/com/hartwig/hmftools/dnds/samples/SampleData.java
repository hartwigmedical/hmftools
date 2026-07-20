package com.hartwig.hmftools.dnds.samples;

import com.hartwig.hmftools.dnds.SampleMutationalLoad;
import com.hartwig.hmftools.dnds.SomaticVariant;

import java.util.List;

public class SampleData {
    public final List<SomaticVariant> Variants;
    public final SampleMutationalLoad MutationalLoad;

    public SampleData(final List<SomaticVariant> variants, final SampleMutationalLoad mutationalLoad) {
        Variants = variants;
        MutationalLoad = mutationalLoad;
    }
}
