package com.hartwig.hmftools.pavereverse.serve;

import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ReversePave;
import com.hartwig.hmftools.pavereverse.variants.ProteinVariant;

class DifferenceWithTransvar
{
    public final BaseSequenceVariants mVariant;
    public final ProteinAnnotationCollator mCollator;
    public final String mType;

    DifferenceWithTransvar(VariantStatus variantStatus, ReversePave baseSequenceVariantsCalculator)
    {
        this.mVariant = variantStatus.Variant;
        this.mCollator = variantStatus.Collator;
        ProteinVariant transvalVariant =
                baseSequenceVariantsCalculator.variationParser()
                        .parseGeneVariants(mCollator.Gene, mCollator.Annotation)
                        .iterator()
                        .next();
        mType = transvalVariant.getClass().getSimpleName();
        System.out.println(mType + " difference for: " + mCollator.Gene + " " + mCollator.Annotation
                + ", calculated: " + mVariant.transcriptName() + " canonical: " + mVariant.Transcript.IsCanonical
                + " " + mVariant.changes()
                + " from transvar: " + mCollator.ChangeSequences);
    }
}
