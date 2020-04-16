package com.hartwig.hmftools.isofox.novel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ChimericRead
{
    public final String Id;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;
    public final Cigar Cigar;
    public final String MateChromosome;
    public final int InsertSize;

    private final Map<RegionMatchType,List<TransExonRef>> mTransExonRefs;

    private final int mFlags;
    // public final List<long[]> mMappedCoords;

    public static ChimericRead from(final SAMRecord record)
    {
        return new ChimericRead(
                record.getReadName(), record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getCigar(), record.getMateReferenceName(), record.getFlags(), record.getInferredInsertSize());
    }

    public ChimericRead(
            final String id, final String chromosome, long posStart, long posEnd, @NotNull final Cigar cigar,
            final String mateChromosome, int flags, int insertSize)
    {
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        Cigar = cigar;
        mFlags = flags;
        MateChromosome = mateChromosome;
        InsertSize = insertSize;

        mTransExonRefs = Maps.newHashMap();
    }

    public final Map<RegionMatchType,List<TransExonRef>> getTransExonRefs() { return mTransExonRefs; }

    public boolean isNegStrand() { return (mFlags & SAMFlag.MATE_REVERSE_STRAND.intValue()) != 0; }
    public boolean isProperPair() { return (mFlags & SAMFlag.PROPER_PAIR.intValue()) != 0; }
    public boolean isDuplicate() { return (mFlags & SAMFlag.DUPLICATE_READ.intValue()) != 0; }
    public boolean isSupplementaryAlignment() { return (mFlags & SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()) != 0; }
    public boolean isSecondaryAlignment() { return (mFlags & SAMFlag.NOT_PRIMARY_ALIGNMENT.intValue()) != 0; }

}
