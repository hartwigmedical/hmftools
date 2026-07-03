package com.hartwig.hmftools.tars.liftback.supplementary;

// Classifies donor/acceptor 2-base flanks: GT-AG canonical (~99%), GC-AG/AT-AC semi-canonical.
// Strand is unknown at scan time so reverse-complement motifs (CT-AC, CT-GC, GT-AT) are accepted.
public final class SpliceMotif
{
    private SpliceMotif() { }

    public static Tier classify(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2)
        {
            return Tier.NONE;
        }
        if(acceptorBases == null || acceptorBases.length != 2)
        {
            return Tier.NONE;
        }

        char d0 = upper(donorBases[0]);
        char d1 = upper(donorBases[1]);
        char a0 = upper(acceptorBases[0]);
        char a1 = upper(acceptorBases[1]);

        if((d0 == 'G' && d1 == 'T' && a0 == 'A' && a1 == 'G')
                || (d0 == 'C' && d1 == 'T' && a0 == 'A' && a1 == 'C'))
            return Tier.CANONICAL;

        if((d0 == 'G' && d1 == 'C' && a0 == 'A' && a1 == 'G')
                || (d0 == 'C' && d1 == 'T' && a0 == 'G' && a1 == 'C'))
            return Tier.SEMI_CANONICAL;

        if((d0 == 'A' && d1 == 'T' && a0 == 'A' && a1 == 'C')
                || (d0 == 'G' && d1 == 'T' && a0 == 'A' && a1 == 'T'))
            return Tier.SEMI_CANONICAL;

        return Tier.NONE;
    }

    private static char upper(final byte b)
    {
        return (char) (b & ~0x20);
    }
}
