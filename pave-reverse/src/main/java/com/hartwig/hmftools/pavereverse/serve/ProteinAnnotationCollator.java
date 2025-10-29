package com.hartwig.hmftools.pavereverse.serve;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

class ProteinAnnotationCollator
{
    final String Chromosome;
    final String Gene;
    final String Annotation;
    Set<BaseSequenceChange> ChangeSequences = new HashSet<>();

    public ProteinAnnotationCollator(ServeItem serveItem)
    {
        Chromosome = serveItem.Chromosome;
        Gene = serveItem.Gene;
        Annotation = serveItem.Annotation;
    }

    public void addHotspot(ServeItem serveItem)
    {
        Preconditions.checkArgument(serveItem.Chromosome.equals(Chromosome));
        Preconditions.checkArgument(serveItem.Gene.equals(Gene));
        Preconditions.checkArgument(serveItem.Annotation.equals(Annotation));
        ChangeSequences.add(new BaseSequenceChange(serveItem.Ref, serveItem.Alt, Chromosome, serveItem.Position));
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ProteinAnnotationCollator that = (ProteinAnnotationCollator) o;
        return Objects.equals(Chromosome, that.Chromosome) && Objects.equals(Gene, that.Gene)
                && Objects.equals(Annotation, that.Annotation);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Gene, Annotation);
    }

    @Override
    public String toString()
    {
        return "ProteinAnnotationCollator{" +
                "mChromosome='" + Chromosome + '\'' +
                ", mGene='" + Gene + '\'' +
                ", mAnnotation='" + Annotation + '\'' +
                '}';
    }
}
