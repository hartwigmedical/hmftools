package com.hartwig.hmftools.pavereverse.serve;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;

class VariantStatus
{
    static VariantStatus withProcessingException(ProteinAnnotationCollator collator, Throwable e)
    {
        VariantStatus status = new VariantStatus((collator));
        status.ProcessingError = e;
        return status;
    }

    final ProteinAnnotationCollator Collator;
    Exception ParsingError = null;
    Throwable ProcessingError = null;
    BaseSequenceVariants Variant;

    VariantStatus(ProteinAnnotationCollator collator, BaseSequenceVariants variant)
    {
        Collator = collator;
        Variant = variant;
    }

    VariantStatus(ProteinAnnotationCollator collator)
    {
        Collator = collator;
    }

    boolean usesNonCanonicalTranscript()
    {
        return !Variant.Transcript.IsCanonical;
    }

    boolean parsedOk()
    {
        return ParsingError == null;
    }

    boolean hasProcessingError()
    {
        return ProcessingError != null;
    }

    boolean haveSameChanges()
    {
        var common = Sets.intersection(Collator.ChangeSequences, Variant.changes());
        if(!common.isEmpty())
        {
            return true;
        }
        if(Collator.Annotation.contains("ins"))
        {
            if(Collator.ChangeSequences.size() == 1 && Variant.changes().size() == 1)
            {
                return changesSameModuloAlt();
            }
        }

        if(Collator.Annotation.endsWith("dup"))
        {
            if(Collator.ChangeSequences.size() == 1 && Variant.changes().size() == 1)
            {
                return changesSameModuloAlt();
            }

            return Collator.ChangeSequences.containsAll(Variant.changes());
        }
        if(Collator.Annotation.endsWith("fs"))
        {
            if(Collator.ChangeSequences.containsAll(Variant.changes()))
            {
                return true;
            }
        }
        return false;
    }

    boolean changesSameModuloAlt()
    {
        if(Variant.changes().size() != 1)
        {
            return false;
        }
        if(Collator.ChangeSequences.size() != 1)
        {
            return false;
        }
        BaseSequenceChange calculatedHS = Variant.changes().iterator().next();
        BaseSequenceChange givenHS = Collator.ChangeSequences.iterator().next();
        if(calculatedHS.Position != givenHS.Position)
        {
            return false;
        }
        if(!calculatedHS.Ref.equals(givenHS.Ref))
        {
            return false;
        }
        return givenHS.Chromosome.equals(calculatedHS.Chromosome);
    }
}
