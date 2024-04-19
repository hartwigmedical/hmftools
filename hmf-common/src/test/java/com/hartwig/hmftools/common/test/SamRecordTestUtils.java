package com.hartwig.hmftools.common.test;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome.MT_LENGTH;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public final class SamRecordTestUtils
{
    public static final int DEFAULT_BASE_QUAL = 37;
    public static final int DEFAULT_MAP_QUAL = 60;

    public static SAMSequenceDictionary SAM_DICTIONARY_V37;

    static
    {
        SAM_DICTIONARY_V37 = new SAMSequenceDictionary();

        RefGenomeCoordinates v37Coords = RefGenomeCoordinates.COORDS_37;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            SAM_DICTIONARY_V37.addSequence(new SAMSequenceRecord(chromosome.toString(), v37Coords.Lengths.get(chromosome)));
        }

        SAM_DICTIONARY_V37.addSequence(new SAMSequenceRecord(MitochondrialChromosome.MT.toString(), MT_LENGTH));
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment,
            boolean mateReversed, final String mateCigar)
    {
        SAMRecord record = createSamRecord(readId, chrStr, readStart, readBases, cigar, mateChr, mateStart, isReversed, isSupplementary, suppAlignment);

        if(mateReversed)
            record.setMateNegativeStrandFlag(true);

        if(mateCigar != null)
            record.setAttribute(MATE_CIGAR_ATTRIBUTE, mateCigar);

        return record;
    }

    public static SAMRecord createSamRecordUnpaired(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, boolean isReversed,
            boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecord record = createSamRecord(
                readId, chrStr, readStart, readBases, cigar, NO_CHROMOSOME_NAME, NO_POSITION, isReversed, isSupplementary, suppAlignment);
        record.setReadPairedFlag(false);
        return record;
    }

    public static SAMRecord cloneSamRecord(final SAMRecord record, final String newReadId)
    {
        try
        {
            SAMRecord newRecord = (SAMRecord)record.clone();
            newRecord.setReadName(newReadId);
            return newRecord;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.getHeader().setSequenceDictionary(SAM_DICTIONARY_V37);
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        int chrIndex = !chrStr.equals(NO_CHROMOSOME_NAME) ? chromosomeOrdinal(chrStr) : NO_CHROMOSOME_INDEX;

        SAMRecord record = recordBuilder.addFrag(
                readId, chrIndex, readStart, isReversed, false,
                cigar, readBases, DEFAULT_BASE_QUAL, false);

        record.setReadBases(readBases.getBytes());

        byte[] qualities = buildDefaultBaseQuals(readBases.length());

        for(int i = 0; i < readBases.length(); ++i)
        {
            qualities[i] = DEFAULT_BASE_QUAL;
        }

        record.setBaseQualities(qualities);
        record.setReferenceName(chrStr);
        record.setReferenceIndex(chrIndex); // need to override since no header is present
        record.setReadPairedFlag(true);

        if(chrIndex == NO_CHROMOSOME_INDEX)
            record.setReadUnmappedFlag(true);

        if(!mateChr.isEmpty() && !mateChr.equals(NO_CHROMOSOME_NAME))
        {
            record.setMateReferenceName(mateChr);
            record.setMateAlignmentStart(mateStart);
            record.setMateReferenceIndex(chromosomeOrdinal(mateChr));
            record.setProperPairFlag(true);
        }
        else
        {
            record.setMateUnmappedFlag(true);
            record.setMateReferenceName(NO_CHROMOSOME_NAME);
            record.setMateAlignmentStart(readStart);
            record.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
            record.setProperPairFlag(false);
        }

        // to be correct this should match the cigar element count
        record.setFirstOfPairFlag(true);

        record.setSupplementaryAlignmentFlag(isSupplementary);

        if(suppAlignment != null)
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppAlignment.asSamTag());

        if(chrStr.equals(mateChr))
            record.setInferredInsertSize(abs(readStart - mateStart));
        else
            record.setInferredInsertSize(0);

        return record;
    }

    public static int chromosomeOrdinal(final String chromosome)
    {
        if(HumanChromosome.contains(chromosome))
            return HumanChromosome.fromString(chromosome).ordinal();

        if(MitochondrialChromosome.contains(chromosome))
            return HumanChromosome.values().length;

        else return HumanChromosome.values().length + 1;
    }

    public static byte[] buildDefaultBaseQuals(int length) { return buildBaseQuals(length, DEFAULT_BASE_QUAL); }

    public static byte[] buildBaseQuals(int length, int qualValue)
    {
        final byte[] baseQuals = new byte[length];

        for(int i = 0; i < length; ++i)
        {
            baseQuals[i] = (byte)qualValue;
        }

        return baseQuals;
    }

    public static int setReadFlag(int flags, final SAMFlag flag)
    {
        flags |= flag.intValue();
        return flags;
    }
}
