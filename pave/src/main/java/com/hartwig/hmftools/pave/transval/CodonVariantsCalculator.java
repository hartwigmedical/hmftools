package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.Checks.isValidAminoAcidName;

import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.AminoAcids;

import org.jetbrains.annotations.NotNull;

public class CodonVariantsCalculator
{
    @NotNull public final String Reference;
    @NotNull public final String Variant;

    public CodonVariantsCalculator(@NotNull final String reference, @NotNull final String variant)
    {
        Preconditions.checkArgument(isValidAminoAcidName(reference));
        Preconditions.checkArgument(isValidAminoAcidName(variant));
        Reference = AminoAcids.forceSingleLetterProteinAnnotation(reference);
        Variant = AminoAcids.forceSingleLetterProteinAnnotation(variant);
    }

    @NotNull
    public SortedSet<CodonChange> possibleCodonsForVariant(@NotNull final String referenceCodon)
    {
        Preconditions.checkArgument(Objects.equals(AminoAcids.findAminoAcidForCodon(referenceCodon), Reference));
        SortedSet<CodonChange> result = new TreeSet<>();
        AminoAcids.AMINO_ACID_TO_CODON_MAP.get(Variant).forEach(variant -> result.add(new CodonChange(referenceCodon, variant)));
        return result;
    }
}
