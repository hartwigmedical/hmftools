package com.hartwig.hmftools.pavereverse.serve;

import java.util.Objects;

import com.google.common.base.Preconditions;

class ServeItem
{
    public final String Gene;
    public final String Annotation;
    public final String Chromosome;
    public final String Ref;
    public final String Alt;
    public final int Position;

    ServeItem(String gene, String annotation, String chromosome, String ref, String alt, final int position)
    {
        Preconditions.checkArgument(!gene.isEmpty());
        Preconditions.checkArgument(!annotation.isEmpty());
        Gene = gene;
        Annotation = annotation;
        Chromosome = chromosome;
        Ref = ref;
        Alt = alt;
        Position = position;
    }

    @Override
    public String toString()
    {
        return "ServeItem{" +
                "Gene='" + Gene + '\'' +
                ", Annotation='" + Annotation + '\'' +
                ", Chromosome='" + Chromosome + '\'' +
                ", Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", Position=" + Position +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ServeItem serveItem = (ServeItem) o;
        return Position == serveItem.Position && Objects.equals(Gene, serveItem.Gene)
                && Objects.equals(Annotation, serveItem.Annotation) && Objects.equals(Chromosome, serveItem.Chromosome)
                && Objects.equals(Ref, serveItem.Ref) && Objects.equals(Alt, serveItem.Alt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Gene, Annotation, Chromosome, Ref, Alt, Position);
    }
}
