package com.hartwig.hmftools.qsee.prep.bqr;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.BqrRecord;

public class ExtendedBqrRecord
{
    public final ConsensusType ReadType;
    public final String StandardMutation;
    public final String StandardTrinucContext;
    public final byte OriginalQuality;
    public final double RecalibratedQuality;
    public final int Count;

    public ExtendedBqrRecord(ConsensusType readType, char refBase, char altBase, String trinucContext, int count, byte originalQual,
            double recalibratedQual)
    {
        ReadType = readType;
        OriginalQuality = originalQual;
        RecalibratedQuality = recalibratedQual;
        Count = count;

        char standardRefBase = refBase;
        char standardAltBase = altBase;
        if(refBase == 'G' || refBase == 'A')
        {
            standardRefBase = swapDnaBase(refBase);
            standardAltBase = swapDnaBase(altBase);
        }

        String standardMutation = String.format("%s>%s", standardRefBase, standardAltBase);

        String standardTrinucContext = (refBase == standardRefBase)
                ? trinucContext
                : reverseComplementBases(trinucContext);

        StandardMutation = standardMutation;
        StandardTrinucContext = standardTrinucContext;
    }

    public static ExtendedBqrRecord from(BqrRecord bqrRecord)
    {
        return new ExtendedBqrRecord(bqrRecord.Key.ReadType,
                (char) bqrRecord.Key.Ref,
                (char) bqrRecord.Key.Alt,
                new String(bqrRecord.Key.TrinucleotideContext),
                bqrRecord.Count,
                bqrRecord.Key.Quality,
                bqrRecord.RecalibratedQuality);
    }
}
