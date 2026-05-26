package com.hartwig.hmftools.redux.splice.rescue;

// Classifies a candidate intron's 2-base donor/acceptor flanks. ~99% of human introns are GT-AG
// (canonical), ~0.6% are GC-AG or AT-AC (semi-canonical). Read strand is unknown at scan time so
// we accept the reverse-strand complement (CT-AC, CT-GC, GT-AT) on equal terms.
public final class SpliceMotif
{
    public static final int TIER_NONE = 0;
    public static final int TIER_SEMI_CANONICAL = 1;
    public static final int TIER_CANONICAL = 2;
    public static final int TIER_ANNOTATED = 3;

    private SpliceMotif() {}

    // donor bases are at [intronStart, intronStart + 1]; acceptor bases at [intronEnd - 1, intronEnd].
    // null or wrong-length inputs return TIER_NONE.
    public static int classify(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2)
            return TIER_NONE;
        if(acceptorBases == null || acceptorBases.length != 2)
            return TIER_NONE;

        final char d0 = upper(donorBases[0]);
        final char d1 = upper(donorBases[1]);
        final char a0 = upper(acceptorBases[0]);
        final char a1 = upper(acceptorBases[1]);

        if((d0 == 'G' && d1 == 'T' && a0 == 'A' && a1 == 'G')
                || (d0 == 'C' && d1 == 'T' && a0 == 'A' && a1 == 'C'))
            return TIER_CANONICAL;

        if((d0 == 'G' && d1 == 'C' && a0 == 'A' && a1 == 'G')
                || (d0 == 'C' && d1 == 'T' && a0 == 'G' && a1 == 'C'))
            return TIER_SEMI_CANONICAL;

        if((d0 == 'A' && d1 == 'T' && a0 == 'A' && a1 == 'C')
                || (d0 == 'G' && d1 == 'T' && a0 == 'A' && a1 == 'T'))
            return TIER_SEMI_CANONICAL;

        return TIER_NONE;
    }

    private static char upper(final byte b)
    {
        return (char)(b & ~0x20);
    }
}
