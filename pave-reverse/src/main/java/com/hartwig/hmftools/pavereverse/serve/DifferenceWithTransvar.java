package com.hartwig.hmftools.pavereverse.serve;

import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ProteinVariant;
import com.hartwig.hmftools.pavereverse.ReversePave;

class DifferenceWithTransvar
{
    public final BaseSequenceVariants variant;
    public final ProteinAnnotationCollator collator;
    public final String type;

    DifferenceWithTransvar(final VariantStatus variantStatus, ReversePave baseSequenceVariantsCalculator)
    {
        this.variant = variantStatus.variant;
        this.collator = variantStatus.collator;
        ProteinVariant transvalVariant =
                baseSequenceVariantsCalculator.variationParser().parseGeneVariants(collator.mGene, collator.mAnnotation).iterator().next();
        type = transvalVariant.getClass().getSimpleName();
        System.out.println(type + " difference for: " + collator.mGene + " " + collator.mAnnotation
                + ", calculated: " + variant.transcriptName() + " canonical: " + variant.mTranscript.IsCanonical
                + " " + variant.changes()
                + " from transvar: " + collator.hotspots);
    }
}
