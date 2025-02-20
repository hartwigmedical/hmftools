package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public class CodonChange implements Comparable<CodonChange>
{
    @NotNull
    public final String ReferenceCodon;
    @NotNull
    public final String AlternateCodon;
    private final int EditDistance;

    public CodonChange(@NotNull final String referenceCodon, @NotNull final String alternateCodon)
    {
        Preconditions.checkArgument(Checks.isCodon(referenceCodon));
        Preconditions.checkArgument(Checks.isCodon(alternateCodon));
        this.ReferenceCodon = referenceCodon;
        this.AlternateCodon = alternateCodon;
        int distance = 0;
        for(int i = 0; i < 3; ++i)
        {
            if(referenceCodon.charAt(i) != alternateCodon.charAt(i))
            {
                distance += 1;
            }
        }
        EditDistance = distance;
    }

    @NotNull
    public CodonChange reverseComplement()
    {
        return new CodonChange(Nucleotides.reverseComplementBases(ReferenceCodon), Nucleotides.reverseComplementBases(AlternateCodon));
    }

    @Override
    public int compareTo(@NotNull final CodonChange o)
    {
        if(!ReferenceCodon.equals(o.ReferenceCodon))
        {
            throw new IllegalArgumentException(String.format("%s != %s", ReferenceCodon, o.ReferenceCodon));
        }
        int byDistance = editDistance() - o.editDistance();
        if(byDistance != 0)
        {
            return byDistance;
        }
        return AlternateCodon.compareTo(o.AlternateCodon);
    }

    public int editDistance()
    {
        return EditDistance;
    }

    public int positionOfFirstDifference()
    {
        for(int i = 0; i < 3; ++i)
        {
            if(ReferenceCodon.charAt(i) != AlternateCodon.charAt(i))
            {
                return i;
            }
        }
        return -1;
    }

    public TransvalHotspot hotspot(GeneData gene, int codonPosition)
    {
        Pair<String, String> refAlt = differenceStrings();
        if(gene.forwardStrand())
        {
            int position = codonPosition + positionOfFirstDifference();
            return new TransvalHotspot(refAlt.getLeft(), refAlt.getRight(), gene.Chromosome, position);
        }
        return new TransvalHotspot(reverseComplementBases(refAlt.getLeft()),
                reverseComplementBases(refAlt.getRight()),
                gene.Chromosome,
                codonPosition - positionOfFirstDifference() + 1);
    }

    public Pair<String,String> differenceStrings()
    {
        StringBuilder refBuilder = new StringBuilder();
        StringBuilder altBuilder = new StringBuilder();
        int differencesAnnotated = 0;
        for(int i = 0; i < 3; ++i)
        {
            if (differencesAnnotated == EditDistance)
            {
                break;
            }
            final char refChar = ReferenceCodon.charAt(i);
            final char altChar = AlternateCodon.charAt(i);
            boolean differentHere = refChar != altChar;
            if (differentHere || differencesAnnotated > 0)
            {
                refBuilder.append(refChar);
                altBuilder.append(altChar);
            }
            if (differentHere)
            {
                differencesAnnotated++;
            }
        }
        return Pair.of(refBuilder.toString(), altBuilder.toString());
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CodonChange that = (CodonChange) o;
        return Objects.equals(ReferenceCodon, that.ReferenceCodon) && Objects.equals(AlternateCodon, that.AlternateCodon);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(ReferenceCodon, AlternateCodon);
    }

    @Override
    public String toString()
    {
        return "CodonChange{" +
                "referenceCodon='" + ReferenceCodon + '\'' +
                ", alternateCodon='" + AlternateCodon + '\'' +
                '}';
    }
}
