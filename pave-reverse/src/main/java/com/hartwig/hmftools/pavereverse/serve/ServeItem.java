package com.hartwig.hmftools.pavereverse.serve;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

class ServeItem
{
    @NotNull
    public final String Gene;

    @NotNull
    public final String Annotation;

    @NotNull
    public final String Chromosome;

    @NotNull
    public final String Ref;

    @NotNull
    public final String Alt;

    public final int Position;

    ServeItem(@NotNull final String gene, @NotNull final String annotation, @NotNull final String chromosome, @NotNull final String ref,
            @NotNull final String alt, final int position)
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
