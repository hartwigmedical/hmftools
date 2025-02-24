package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Lists;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.codon.Codons;

import org.jetbrains.annotations.NotNull;

class AminoAcidSequence
{
    @NotNull
    private final List<AminoAcid> aminoAcids;

    @NotNull
    public static AminoAcidSequence parse(@NotNull final String sequence)
    {
        Pattern variantListPattern = Pattern.compile("([A-Z][a-z]*)");
        final Matcher variantListMatcher = variantListPattern.matcher(sequence);
        List<AminoAcid> aminoAcids = Lists.newArrayList();
        while(variantListMatcher.find())
        {
            String aa = variantListMatcher.group();
            aminoAcids.add(new AminoAcid(aa));
        }

        return new AminoAcidSequence(aminoAcids);
    }

    public static AminoAcidSequence fromNucleotides(@NotNull final String nucleotides)
    {
        Preconditions.checkArgument(!nucleotides.isEmpty(), "Non-empty nucleotide sequence required");
        Preconditions.checkArgument(nucleotides.length() % 3 == 0, "Length of sequence must be a multiple of 3");
        List<AminoAcid> aminoAcids = Lists.newArrayList();
        for (int i = 0; i < nucleotides.length(); i += 3)
        {
            String codon = nucleotides.substring(i, i + 3);
            String acid = AminoAcids.findAminoAcidForCodon(codon);
            if(acid != null)
            {
                aminoAcids.add(new AminoAcid(acid));
            }
            else if(Codons.isStopCodon(codon))
            {
                aminoAcids.add(new AminoAcid("X"));
                break;
            }
            else
            {
                throw new IllegalArgumentException("Unknown codon: " + codon);
            }
        }
        return new AminoAcidSequence(aminoAcids);
    }

    public static AminoAcidSequence empty()
    {
        return new AminoAcidSequence(Lists.newArrayList());
    }

    AminoAcidSequence(@NotNull final List<AminoAcid> aminoAcids)
    {
        this.aminoAcids = List.copyOf(aminoAcids);
    }

    public String symbolAt(int index)
    {
        return get(index).symbol;
    }

    public AminoAcid get(int index)
    {
        Preconditions.checkArgument(index >= 0 && index < aminoAcids.size(), "Index out of bounds");
        return aminoAcids.get(index);
    }

    public String sequence()
    {
        return aminoAcids.stream().map(aa -> aa.symbol).collect(Collectors.joining());
    }

    public int length()
    {
        return aminoAcids.size();
    }

    @NotNull
    public AminoAcidSequence subsequenceUpToInclusive(int i)
    {
        Preconditions.checkArgument(i >= 0 && i <= aminoAcids.size(), "Index out of bounds");
        return new AminoAcidSequence(aminoAcids.subList(0, i));
    }

    @NotNull
    public AminoAcidSequence subsequenceAfterExclusive(int i)
    {
        Preconditions.checkArgument(i >= 0 && i < aminoAcids.size(), "Index out of bounds");
        return new AminoAcidSequence(aminoAcids.subList(i + 1, aminoAcids.size()));
    }

    @NotNull
    public AminoAcidSequence append(@NotNull AminoAcidSequence other)
    {
        List<AminoAcid> newAminoAcids = Lists.newArrayList(aminoAcids);
        newAminoAcids.addAll(other.aminoAcids);
        return new AminoAcidSequence(newAminoAcids);
    }

    @NotNull
    public AminoAcidSequence duplicateRange(int startInclusive, int endExclusive)
    {
        List<AminoAcid> newAminoAcids = Lists.newArrayList(aminoAcids.subList(0, startInclusive));
        newAminoAcids.addAll(aminoAcids.subList(startInclusive, endExclusive));
        newAminoAcids.addAll(aminoAcids.subList(startInclusive, endExclusive));
        newAminoAcids.addAll(aminoAcids.subList(endExclusive, aminoAcids.size()));
        return new AminoAcidSequence(newAminoAcids);
    }

    @NotNull
    public AminoAcidSequence deleteRange(int startInclusive, int endExclusive)
    {
        List<AminoAcid> newAminoAcids = Lists.newArrayList(aminoAcids.subList(0, startInclusive));
        newAminoAcids.addAll(aminoAcids.subList(endExclusive, aminoAcids.size()));
        return new AminoAcidSequence(newAminoAcids);
    }

    @NotNull
    public AminoAcidSequence insert(int positionAfterWhichToInsert, AminoAcidSequence toInsert)
    {
        List<AminoAcid> newAminoAcids = Lists.newArrayList(aminoAcids.subList(0, positionAfterWhichToInsert));
        newAminoAcids.addAll(toInsert.aminoAcids);
        newAminoAcids.addAll(aminoAcids.subList(positionAfterWhichToInsert, aminoAcids.size()));
        return new AminoAcidSequence(newAminoAcids);
    }

    @NotNull
    public AminoAcidSequence replace(int position, AminoAcid replacement)
    {
        Preconditions.checkArgument(position >= 1 && position <= aminoAcids.size(), "Index out of bounds");
        List<AminoAcid> newAminoAcids = Lists.newArrayList(aminoAcids.subList(0, position - 1));
        newAminoAcids.add(replacement);
        newAminoAcids.addAll(aminoAcids.subList(position, aminoAcids.size()));
        return new AminoAcidSequence(newAminoAcids);
    }

    @Override
    public String toString()
    {
        return sequence();
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final AminoAcidSequence that = (AminoAcidSequence) o;
        return Objects.equals(aminoAcids, that.aminoAcids);
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(aminoAcids);
    }
}
