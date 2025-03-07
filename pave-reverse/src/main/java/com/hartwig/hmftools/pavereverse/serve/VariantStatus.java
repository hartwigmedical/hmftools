package com.hartwig.hmftools.pavereverse.serve;

import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;

import org.jetbrains.annotations.NotNull;

class VariantStatus
{
    static VariantStatus withProcessingException(@NotNull final ProteinAnnotationCollator collator, final @NotNull Throwable e)
    {
        VariantStatus status = new VariantStatus((collator));
        status.mProcessingError = e;
        return status;
    }

    @NotNull
    final ProteinAnnotationCollator collator;

    Exception mParseException = null;

    Throwable mProcessingError = null;

    BaseSequenceVariants variant;

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, BaseSequenceVariants variant)
    {
        this.collator = collator;
        this.variant = variant;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator)
    {
        this.collator = collator;
    }

    boolean usesNonCanonicalTranscript()
    {
        return !variant.mTranscript.IsCanonical;
    }

    boolean parsedOk()
    {
        return mParseException == null;
    }

    boolean hasProcessingError()
    {
        return mProcessingError != null;
    }

    boolean hotspotsSame()
    {
        var common = Sets.intersection(collator.hotspots, variant.changes());
        if(!common.isEmpty())
        {
            return true;
        }
        if(collator.mAnnotation.contains("ins"))
        {
            if(collator.hotspots.size() == 1 && variant.changes().size() == 1)
            {
                return hotspotsSameModuloAlt();
            }
        }

        if(collator.mAnnotation.endsWith("dup"))
        {
            if(collator.hotspots.size() == 1 && variant.changes().size() == 1)
            {
                return hotspotsSameModuloAlt();
            }

            return collator.hotspots.containsAll(variant.changes());
        }
        if(collator.mAnnotation.endsWith("fs"))
        {
            if(collator.hotspots.containsAll(variant.changes()))
            {
                return true;
            }
        }
        return false;
    }

    boolean hotspotsSameModuloAlt()
    {
        if(variant.changes().size() != 1)
        {
            return false;
        }
        if(collator.hotspots.size() != 1)
        {
            return false;
        }
        BaseSequenceChange calculatedHS = variant.changes().iterator().next();
        BaseSequenceChange givenHS = collator.hotspots.iterator().next();
        if(calculatedHS.mPosition != givenHS.mPosition)
        {
            return false;
        }
        if(!calculatedHS.Ref.equals(givenHS.Ref))
        {
            return false;
        }
        return givenHS.mChromosome.equals(calculatedHS.mChromosome);
    }
}
