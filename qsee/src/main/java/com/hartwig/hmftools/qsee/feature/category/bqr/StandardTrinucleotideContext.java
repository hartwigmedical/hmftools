package com.hartwig.hmftools.qsee.feature.category.bqr;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import com.hartwig.hmftools.common.genome.region.Orientation;

public record StandardTrinucleotideContext(
        char Ref,
        char Alt,
        String TrinucleotideContext,
        Orientation ReadOrientation
){
    /*
    The forward (aka standard) context has the ref base as a C or T

    Example target forward mutation context: GTA_T>G

    1) observed, fw orientation: 5' ...G[T>G]A... 3'
       complementary strand:     3' ...C[A>C]T... 5'

    2) observed, rv orientation: 5' ...T[A>C]C... 3'  (complementary strand reversed)
     */

    public static StandardTrinucleotideContext fromNucleotides(char refBase, char altBase, String trinucleotideContext)
    {

        char standardRefBase;
        char standardAltBase;

        boolean refBaseIsStandard = refBase == 'C' || refBase == 'T';

        if(refBaseIsStandard)
        {
            standardRefBase = refBase;
            standardAltBase = altBase;
        }
        else
        {
            standardRefBase = swapDnaBase(refBase);
            standardAltBase = swapDnaBase(altBase);
        }

        Orientation readOrientation = (refBase == standardRefBase)
                ? Orientation.FORWARD
                : Orientation.REVERSE;

        String standardTrinucleotideContext = readOrientation.isForward()
                ? trinucleotideContext
                : reverseComplementBases(trinucleotideContext);

        return new StandardTrinucleotideContext(
                standardRefBase,
                standardAltBase,
                standardTrinucleotideContext,
                readOrientation
        );
    }

    public String mutationContextString()
    {
        return String.format("%s[%s>%s]%s %s",
                TrinucleotideContext.charAt(0), Ref, Alt, TrinucleotideContext.charAt(2), ReadOrientation.asChar());
    }

    public String mutationString()
    {
        String mutation = String.format("%s>%s", Ref, Alt);

        if(Ref == 'C' && Alt == 'T')
        {
            boolean isCpG = TrinucleotideContext.charAt(2) == 'G';
            mutation += isCpG ? " CpG" : " other";
        }

        return mutation;
    }
}
