package com.hartwig.hmftools.qsee.feature.category.bqr;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.BqrKey;
import com.hartwig.hmftools.common.redux.BqrRecord;

public record StandardBqrRecord(
        StandardTrinucleotideContext StandardContext,
        byte OriginalQuality,
        ConsensusType ReadType,
        int Count,
        double RecalibratedQuality
){
    public static StandardBqrRecord from(BqrRecord bqrRecord)
    {
        BqrKey key = bqrRecord.Key;

        StandardTrinucleotideContext standardContext = StandardTrinucleotideContext.fromNucleotides(
                (char) key.Ref,
                (char) key.Alt,
                new String(key.TrinucleotideContext)
        );

        return new StandardBqrRecord(standardContext, key.Quality, key.ReadType, bqrRecord.Count, bqrRecord.RecalibratedQuality);
    }
}
