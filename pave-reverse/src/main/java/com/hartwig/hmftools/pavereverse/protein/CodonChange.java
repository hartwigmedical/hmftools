package com.hartwig.hmftools.pavereverse.protein;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.pavereverse.util.Checks;
import com.hartwig.hmftools.pavereverse.util.PRUtils;

import org.apache.commons.lang3.tuple.Pair;

public class CodonChange implements Comparable<CodonChange>
{
    public final String ReferenceCodon;
    public final String AlternateCodon;
    private final int EditDistance;

    public CodonChange(String referenceCodon, String alternateCodon)
    {
        Preconditions.checkArgument(Checks.isCodon(referenceCodon));
        Preconditions.checkArgument(Checks.isCodon(alternateCodon));
        ReferenceCodon = referenceCodon;
        AlternateCodon = alternateCodon;
        EditDistance = PRUtils.substitutionDistance(referenceCodon, alternateCodon);
    }

    public CodonChange reverseComplement()
    {
        return new CodonChange(Nucleotides.reverseComplementBases(ReferenceCodon), Nucleotides.reverseComplementBases(AlternateCodon));
    }

    @Override
    public int compareTo(CodonChange o)
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

    public Pair<String, String> differenceStrings()
    {
        StringBuilder refBuilder = new StringBuilder();
        StringBuilder altBuilder = new StringBuilder();
        int differencesAnnotated = 0;
        for(int i = 0; i < 3; ++i)
        {
            if(differencesAnnotated == EditDistance)
            {
                break;
            }
            final char refChar = ReferenceCodon.charAt(i);
            final char altChar = AlternateCodon.charAt(i);
            boolean differentHere = refChar != altChar;
            if(differentHere || differencesAnnotated > 0)
            {
                refBuilder.append(refChar);
                altBuilder.append(altChar);
            }
            if(differentHere)
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
                "ReferenceCodon='" + ReferenceCodon + '\'' +
                ", AlternateCodon='" + AlternateCodon + '\'' +
                '}';
    }
}
