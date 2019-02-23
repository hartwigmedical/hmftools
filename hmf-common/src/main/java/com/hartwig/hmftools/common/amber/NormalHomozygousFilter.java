package com.hartwig.hmftools.common.amber;

import java.util.function.Predicate;

public class NormalHomozygousFilter implements Predicate<NormalBAF> {
    @Override
    public boolean test(final NormalBAF bafEvidence) {
        return bafEvidence.isValid(1) && bafEvidence.altSupport() == 0;
    }

}
