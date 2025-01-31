package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Lists;

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

    AminoAcidSequence(@NotNull final List<AminoAcid> aminoAcids)
    {
        this.aminoAcids = List.copyOf(aminoAcids);
    }

    public String sequence()
    {
        return aminoAcids.stream().map(aa -> aa.symbol).collect(Collectors.joining());
    }
}
