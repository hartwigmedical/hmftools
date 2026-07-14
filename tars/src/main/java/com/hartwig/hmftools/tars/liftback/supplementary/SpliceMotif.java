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

        char donor0 = upper(donorBases[0]);
        char donor1 = upper(donorBases[1]);
        char acceptor0 = upper(acceptorBases[0]);
        char acceptor1 = upper(acceptorBases[1]);

        if((donor0 == 'G' && donor1 == 'T' && acceptor0 == 'A' && acceptor1 == 'G')
                || (donor0 == 'C' && donor1 == 'T' && acceptor0 == 'A' && acceptor1 == 'C'))
            return Tier.CANONICAL;

        if((donor0 == 'G' && donor1 == 'C' && acceptor0 == 'A' && acceptor1 == 'G')
                || (donor0 == 'C' && donor1 == 'T' && acceptor0 == 'G' && acceptor1 == 'C'))
            return Tier.SEMI_CANONICAL;

        if((donor0 == 'A' && donor1 == 'T' && acceptor0 == 'A' && acceptor1 == 'C')
                || (donor0 == 'G' && donor1 == 'T' && acceptor0 == 'A' && acceptor1 == 'T'))
            return Tier.SEMI_CANONICAL;

        return Tier.NONE;
    }

    private static char upper(final byte base)
    {
        return (char) (base & ~0x20);
    }
}
