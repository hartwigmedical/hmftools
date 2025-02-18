package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.Checks.isNucleotideSequence;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

class CodonWithinExons
{
    static CodonWithinExons factory(@NotNull String left, @NotNull BaseSequence body, @NotNull String right)
    {
        Preconditions.checkArgument(left.isEmpty() || right.isEmpty());
        Preconditions.checkArgument(left.isEmpty() || isNucleotideSequence(left));
        Preconditions.checkArgument(right.isEmpty() || isNucleotideSequence(right));
        // todo edge cases
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

    CodonWithinExons(final BaseSequence body)
    {
        Body = body;
    }

    @NotNull
    String fixedPrefix()
    {
        return "";
    }

    @NotNull
    String fixedSuffix()
    {
        return "";
    }

    @NotNull
    String variablePart()
    {
        return Body.Bases;
    }

    int strandLocationOfStartOfVariablePart()
    {
        return Body.Start;
    }

    int strandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= 0);
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2);
        return Body.Start + positionWithinCodonOfChange;
    }

    @NotNull
    public Set<CodonChange> possibleVariantsGiving(AminoAcid aminoAcid) {
        Set<String> possibleCodons = aminoAcid.matchingCodons(fixedPrefix(), fixedSuffix());
        final HashSet<CodonChange> result = new HashSet<>();
        possibleCodons.forEach(altCodon -> {
            result.add(new CodonChange(codon(), altCodon));
        });
        return result;
    }

    @NotNull
    String codon()
    {
        return fixedPrefix() + variablePart() + fixedSuffix();
    }
}

class CodonWithPrefixInPreviousExon extends CodonWithinExons
{
    @NotNull
    private final String Prefix;

    CodonWithPrefixInPreviousExon(@NotNull String prefix, final BaseSequence body)
    {
        super(body);
        this.Prefix = prefix;
    }

    int strandLocationOfChange(int positionWithinCodonOfChange)
    {
        Preconditions.checkArgument(positionWithinCodonOfChange >= Prefix.length());
        Preconditions.checkArgument(positionWithinCodonOfChange <= 2);
        return Body.Start + positionWithinCodonOfChange - Prefix.length();
    }

    @NotNull
    @Override
    String fixedPrefix()
    {
        return Prefix;
    }
}
class CodonWithSuffixInNextExon extends CodonWithinExons
{
    @NotNull
    private final String Suffix;

    CodonWithSuffixInNextExon(final BaseSequence body, @NotNull String suffix)
    {
        super(body);
        this.Suffix = suffix;
    }

    @NotNull
    @Override
    String fixedSuffix()
    {
        return Suffix;
    }
}
