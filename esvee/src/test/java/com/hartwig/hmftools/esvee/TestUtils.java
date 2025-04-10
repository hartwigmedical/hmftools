package com.hartwig.hmftools.esvee;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.setReadFlag;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class TestUtils
{
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_READ_ID_2 = "READ_02";

    public static final String TEST_CIGAR_100  = "100M";

    public static final AssemblyConfig TEST_CONFIG = new AssemblyConfig();

    public static final String REF_BASES_RANDOM_100 = generateRandomBases(100);

    public static final int DEFAULT_MAP_QUAL = 60;
    public static final int DEFAULT_NM = 0;

    public static final ReadIdGenerator READ_ID_GENERATOR = new ReadIdGenerator();

    public static String REF_BASES_200 =
            "AAACCCGGGTTTACGTAACCGGTTACGTAAAAACCCCCGGGGGTTTTTACGTAACCGGTTACGTAAACCCGGGTTTAAACGTTTTTGGGGCCCCAAAAAC"
    //       0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    // index 0         10        20        30        40        50        60        70        80        90
          + "GTGGCCTTAAACGTCCCAAAATTTTGGGGACGTGGGGGCCCCCAAAATTTTACGTCCGGTTAAACGTTTTCCCGGGAAAACGTAACCGGTTGGCCAAGCT";

    public static String REF_BASES_400 =
            "AGTGTCGCGAATGCTGATTCGGAACCTTGATCGTGATGTCGATGGCTAGATCAAATCGTCTAGTGGCTAATGCTCGATCGATATGGCTTGATCGTGAGTC"
    //       0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    // index 0         10        20        30        40        50        60        70        80        90
          + "GGGCTGATGATGTTGATCGTAGTGCTAGTGATCGTAGTTGCCAAATGCTGTTGCCCTGTAGCTGATTGGCTAGCTGTAGCTGATCGATGCAATCCGGTAG"
          + "TGTAGCTGATCGCAGGTCGAACCTGGTGATCGATGTCGATCGACTGATGTAGTAGCTGATCGGATGCATGCGTAGCGATGCTAGCTGATCGATTGGCTAA"
          + "GTCGCTTCCGGTATTTGCGTTCCGGGTTTTTTCCGAGCCTACCCCAGTTGGTTAAAAGGATATTATATATATGGCGGCTATATATGCGGTGTGTGTAACC";

    public static String REF_BASES_600 =
            "ATCATCGAATGGAATGGAATGGAACAGTCAATGAACTCGAATGGAATCATCATTGAATGGAATCGAATGGAATCATCGAGTGGAATCGAATGGAATTATG"
    //       0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    // index 0         10        20        30        40        50        60        70        80        90
          + "ATCAAATGGAATCGAATGTAATCATCATCAAATGGAATCAAAAATAACCATAATATGGTCTTTAAAGAAGCAATCTAGCTAAAATGAAATCATTAATCCA"
    // index 100       110       120       130       140       150       160       170       180       190
          + "ATAGGACTAACATCCTTATGAGAAGACAAGATTAGGACACAGGCATGAACTGGGGGAAGACCATGTAAAGACACTCAAAGAACTGACCGTACCACCCCCT"
    // index 200       210       120       130       140       150       160       170       180       190
          + "TCGAATGAATTGAATGCAATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGTAAAGAAGCATTGAGCTATTTACATAGAAATCTCATTTAACTGTG"
    // index 300       310       120       130       140       150       160       170       180       190
          + "ATATAAATTAACCTTTCTTATCCTGCTTCTAAACAAAGGTAAGGGCCACCCAGTCAATGCTTTGTATTCTTCCAATATTCTTTCCTAGAACTTCTTCAAA"
    // index 400       410       120       130       140       150       160       170       180       190
          + "GGCTCTCATGAAGCACTGGTGAAACTGGAAATCACTGAATTTTACTACCATTTTCTTATCCTGCTTCTAAACAAAGGTAAGTTTCTTATCCTGCTTCTAA";
    // index 500       510       120       130       140       150       160       170       180       190

    public static final MockRefGenome REF_GENOME = new MockRefGenome();

    public static Read createRead(final String readId, int readStart, final String readBases, final String cigar)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, readStart + 1000,
                false, false, null);
        return new Read(record);
    }

    public static Read createRead(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar,
            final String mateChr, int mateStart, boolean mateReversed)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, chromosome, readStart, readBases, cigar, mateChr, mateStart,
                false, false, null);

        if(mateReversed)
            record.setMateNegativeStrandFlag(true);

        record.setAttribute(MATE_CIGAR_ATTRIBUTE, format("%dM", readBases.length()));

        return new Read(record);
    }

    public static void setMateCigar(final Read read, final String mateCigar)
    {
        read.bamRecord().setAttribute(MATE_CIGAR_ATTRIBUTE, mateCigar);
    }

    public static Read createConcordantRead(
            final String readId, int readStart, final String readBases, final String cigar, final int mateStart)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, mateStart,
                false, false, null);

        if(readStart > mateStart)
            record.setReadNegativeStrandFlag(true);
        else
            record.setMateNegativeStrandFlag(true);

        record.setAttribute(MATE_CIGAR_ATTRIBUTE, cigar); // assume the same

        return new Read(record);
    }

    public static Read cloneRead(final Read read, final String newReadId)
    {
        SAMRecord newRecord = cloneSamRecord(read.bamRecord(), newReadId);
        return new Read(newRecord);
    }

    public static String makeCigarString(final String readBases, int scLeft, int scRight)
    {
        StringBuilder sb = new StringBuilder();

        if(scLeft > 0)
            sb.append(format("%dS", scLeft));

        sb.append(format("%dM", readBases.length() - scLeft - scRight));

        if(scRight > 0)
            sb.append(format("%dS", scRight));

        return sb.toString();
    }

    public static void setSecondInPair(final SAMRecord record)
    {
        record.setSecondOfPairFlag(true);
        record.setFirstOfPairFlag(false);
    }


    // SAMRecord convenience methods - reconcile with those in hmf-common
    public static String readIdStr(int readId) { return format("READ_%02d", readId); }

    public static int buildFlags(boolean firstInPair, boolean reversed, boolean supplementary)
    {
        int flags = 0;

        flags = setReadFlag(flags, SAMFlag.READ_PAIRED);
        flags = setReadFlag(flags, SAMFlag.PROPER_PAIR);

        if(reversed)
            flags = setReadFlag(flags, SAMFlag.READ_REVERSE_STRAND);

        if(firstInPair)
            flags = setReadFlag(flags, SAMFlag.FIRST_OF_PAIR);
        else
            flags = setReadFlag(flags, SAMFlag.SECOND_OF_PAIR);

        if(supplementary)
            flags = setReadFlag(flags, SAMFlag.SUPPLEMENTARY_ALIGNMENT);

        return flags;
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar)
    {
        return createSamRecord(
                readId, chromosome, readStart, readBases, cigar,
                buildFlags(true, false, false),
                DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar, int flags)
    {
        return createSamRecord( readId, chromosome, readStart, readBases, cigar, flags, DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String mateChr, int mateStart,
            boolean firstInPair, boolean isSupp, final String suppData)
    {
        String readBases = REF_BASES_200.substring(0, 100);
        SAMRecord record = createSamRecord(
                readId, chromosome, readStart, readBases, "100M",
                buildFlags(firstInPair, false, isSupp),
                DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL);

        record.setMateReferenceName(mateChr);
        record.setMateAlignmentStart(mateStart);

        if(suppData != null && !suppData.isEmpty())
            record.setAttribute(SUPPLEMENTARY_ATTRIBUTE, suppData);

        return record;
    }

    public static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String readBases, final String cigar, int flags,
            int mapQual, int baseQual)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        SAMRecord record = recordBuilder.addFrag(
                readId, 1, readStart, false, false, cigar, readBases, mapQual, false);

        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = (byte)baseQual;

        record.setBaseQualities(qualities);
        record.setReferenceName(chromosome);

        record.setFlags(flags);

        record.setMateReferenceName(chromosome);
        record.setMateAlignmentStart(readStart + 300);
        record.setMateNegativeStrandFlag(true);

        record.setInferredInsertSize(400);

        return record;
    }

    public static int getSupportTypeCount(final JunctionAssembly assembly, final SupportType type)
    {
        return (int)assembly.support().stream().filter(x -> x.type() == type).count();
    }
}