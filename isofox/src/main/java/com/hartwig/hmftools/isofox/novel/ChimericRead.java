package com.hartwig.hmftools.isofox.novel;

import java.util.List;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class ChimericRead
{
    public final String Id;
    public final String Chromosome;
    public final long PosStart;
    public final long PosEnd;

    public final String ReadBases;
    public final int Length; // of bases
    public final Cigar Cigar;

    public final boolean IsDuplicate;
    public final boolean IsFirstOfPair;
    public final boolean IsNegStrand;
    public final String MateChromosome;
    public final boolean MateIsNegStrand;

    public final List<TransExonRef> TransExonData;

    // public final List<long[]> mMappedCoords;

    public static ChimericRead from(final SAMRecord record, final List<TransExonRef> transExonData)
    {
        return new ChimericRead(
                record.getReadName(), record.getReferenceName(), record.getStart(), record.getEnd(),
                record.getReadString(), record.getCigar(), record.getFirstOfPairFlag(),
                record.getReadNegativeStrandFlag(), record.getMateReferenceName(), record.getMateNegativeStrandFlag(),
                record.getDuplicateReadFlag(), transExonData);
    }

    public ChimericRead(
            final String id, final String chromosome, long posStart, long posEnd, final String readBases, @NotNull final Cigar cigar,
            boolean isFirstOfPair, boolean isNegStrand, final String mateChromosome, boolean mateIsNegStrand, boolean isDuplicate,
            final List<TransExonRef> transExonData)
    {
        Id = id;
        Chromosome = chromosome;
        PosStart = posStart;
        PosEnd = posEnd;
        ReadBases = readBases;
        Length = ReadBases.length();
        Cigar = cigar;
        IsFirstOfPair = isFirstOfPair;
        IsNegStrand = isNegStrand;
        MateChromosome = mateChromosome;
        MateIsNegStrand = mateIsNegStrand;
        IsDuplicate = isDuplicate;

        TransExonData = transExonData;

        // mMappedCoords = generateMappedCoords(Cigar, PosStart);
    }

}
