package com.hartwig.hmftools.pavereverse.base;

import static com.hartwig.hmftools.pavereverse.Checks.isNucleotideSequence;

import java.util.HashSet;
import java.util.Set;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.pavereverse.variants.CodonChange;
import com.hartwig.hmftools.pavereverse.aa.AminoAcid;

public class CodonWithinExons
{
    public static CodonWithinExons factory(@NotNull String left, @NotNull BaseSequence body, @NotNull String right)
    {
        Preconditions.checkArgument(left.isEmpty() || right.isEmpty());
        Preconditions.checkArgument(left.isEmpty() || isNucleotideSequence(left));
        Preconditions.checkArgument(right.isEmpty() || isNucleotideSequence(right));
        if(!left.isEmpty())
        {
            return new CodonWithPrefixInPreviousExon(left, body);
        }
        if(!right.isEmpty())
        {
            return new CodonWithSuffixInNextExon(body, right);
        }
        return new CodonWithinExons(body);
    }

    protected final BaseSequence Body;

    public CodonWithinExons(final BaseSequence body)
    {
        Body = body;
    }

    @NotNull
    public String fixedPrefix()
    {
        return "";
    }

    @NotNull
    public String fixedSuffix()
    {
        return "";
    }

    @NotNull
    public String variablePart()
    {
        return Body.Bases;
    }

    public int strandLocationOfStartOfVariablePart()
    {
        return Body.Start;
    }

    public int forwardStrandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= 0);
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2);
        return Body.Start + positionWithinCodonOfChange;
    }

    @NotNull
    public Set<CodonChange> possibleVariantsGivingStop()
    {
        return possibleVariantsGiving(new AminoAcid("X"));
    }

    @NotNull
    public Set<CodonChange> possibleVariantsGiving(AminoAcid aminoAcid)
    {
        Set<String> possibleCodons = aminoAcid.matchingCodons(fixedPrefix(), fixedSuffix());
        final HashSet<CodonChange> result = new HashSet<>();
        possibleCodons.forEach(altCodon -> result.add(new CodonChange(codon(), altCodon)));
        return result;
    }

    @NotNull
    public String codon()
    {
        return fixedPrefix() + variablePart() + fixedSuffix();
    }

    @NotNull
    public CodonWithinExons reverseComplement()
    {
        return new CodonWithinExons(Body.reverseComplement());
    }
}

