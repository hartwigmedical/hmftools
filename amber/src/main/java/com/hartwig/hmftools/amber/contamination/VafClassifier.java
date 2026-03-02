package com.hartwig.hmftools.amber.contamination;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.amber.purity.CanonicalSnvType;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public interface VafClassifier<S, T extends Comparable<T>>
{
    static VafClassifier<VafReading, ChrArm> chrArmClassifier(ChrArmLocator locator)
    {
        return contamination -> locator.map(contamination.chromosome(), contamination.position());
    }

    static VafClassifier<PositionEvidence, ChrArm> chrArmPositionEvidenceClassifier(ChrArmLocator locator)
    {
        return evidence -> locator.map(evidence.chromosome(), evidence.position());
    }

    static VafClassifier<PositionEvidence, CanonicalSnvType> mutationTypeClassifier()
    {
        return reading -> CanonicalSnvType.type(reading.Ref, reading.Alt);
    }

    T classify(S vafReading);
}

