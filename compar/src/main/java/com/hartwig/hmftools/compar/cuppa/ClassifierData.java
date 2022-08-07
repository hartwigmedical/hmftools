package com.hartwig.hmftools.compar.cuppa;

import com.hartwig.hmftools.common.cuppa.ClassifierType;

public class ClassifierData
{
    public final ClassifierType Classifier;
    public final String TopRefCancerType;
    public final double TopRefCancerValue;

    public ClassifierData(final ClassifierType classifier, final String topRefCancerType, final double topRefCancerValue)
    {
        Classifier = classifier;
        TopRefCancerType = topRefCancerType;
        TopRefCancerValue = topRefCancerValue;
    }
}
