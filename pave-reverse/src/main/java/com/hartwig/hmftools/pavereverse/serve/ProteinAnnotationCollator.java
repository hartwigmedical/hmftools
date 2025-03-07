package com.hartwig.hmftools.pavereverse.serve;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.jetbrains.annotations.NotNull;

class ProteinAnnotationCollator
{
    @NotNull
    final String mChromosome;

    @NotNull
    final String mGene;

    @NotNull
    final String mAnnotation;

    @NotNull
    Set<BaseSequenceChange> hotspots = new HashSet<>();

    public ProteinAnnotationCollator(ServeItem serveItem)
    {
        mChromosome = serveItem.Chromosome;
        mGene = serveItem.Gene;
        mAnnotation = serveItem.Annotation;
    }

    public void addHotspot(ServeItem serveItem)
    {
        Preconditions.checkArgument(serveItem.Chromosome.equals(mChromosome));
        Preconditions.checkArgument(serveItem.Gene.equals(mGene));
        Preconditions.checkArgument(serveItem.Annotation.equals(mAnnotation));
        hotspots.add(new BaseSequenceChange(serveItem.Ref, serveItem.Alt, mChromosome, serveItem.Position));
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ProteinAnnotationCollator that = (ProteinAnnotationCollator) o;
        return Objects.equals(mChromosome, that.mChromosome) && Objects.equals(mGene, that.mGene)
                && Objects.equals(mAnnotation, that.mAnnotation);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, mGene, mAnnotation);
    }

    @Override
    public String toString()
    {
        return "ProteinAnnotationCollator{" +
                "mChromosome='" + mChromosome + '\'' +
                ", mGene='" + mGene + '\'' +
                ", mAnnotation='" + mAnnotation + '\'' +
                '}';
    }
}
