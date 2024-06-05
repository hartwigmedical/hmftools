package com.hartwig.hmftools.esvee;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.setReadFlag;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
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

        return new Read(record);
    }

    public static Read createConcordantRead(
            final String readId, int readStart, final String readBases, final String cigar, final int mateStart)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, mateStart,
                false, false, null);

        record.setMateNegativeStrandFlag(true);

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

    public static List<SAMRecord> createJunctionReads(
            final MockRefGenome refGenome, final String readId, int anchorLength,
            final String chrStart, int junctionPosStart, Orientation junctionOrientStart,
            final String chrEnd, int junctionPosEnd, Orientation junctionOrientEnd, int mateStart)
    {
        // creates a junction read, its supplementary and a local mate if the coords are supplied
        int readBaseLength = anchorLength * 2;
        int readStart, readEnd, suppStart, suppEnd;
        String readCigar, suppCigar;
        String basesStart, basesEnd;

        if(junctionOrientStart.isForward())
        {
            readStart = junctionPosStart - anchorLength + 1;
            readEnd = junctionPosStart;
            readCigar = format("%dM%dS", anchorLength, anchorLength);
        }
        else
        {
            readStart = junctionPosStart;
            readEnd = junctionPosStart + anchorLength - 1;
            readCigar = format("%dS%dM", anchorLength, anchorLength);
        }

        basesStart = refGenome.getBaseString(chrStart, readStart, readEnd);

        if(junctionOrientEnd.isForward())
        {
            suppStart = junctionPosEnd - anchorLength + 1;
            suppEnd = junctionPosEnd;
            suppCigar = format("%dM%dS", anchorLength, anchorLength);
        }
        else
        {
            suppStart = junctionPosEnd;
            suppEnd = junctionPosEnd + anchorLength - 1;
            suppCigar = format("%dS%dM", anchorLength, anchorLength);
        }

        basesEnd = refGenome.getBaseString(chrEnd, suppStart, suppEnd);

        String readBases, suppBases;
        boolean isSuppNegStrand = true;

        if(junctionOrientStart != junctionOrientEnd)
        {
            if(junctionOrientStart.isForward())
                readBases = basesStart + basesEnd;
            else
                readBases = basesEnd + basesStart;

            suppBases = readBases;
        }
        else
        {
            isSuppNegStrand = false;

            // keep the first read's bases in the 5' to 3' direction
            if(junctionOrientStart.isForward())
            {
                readBases = basesStart + Nucleotides.reverseComplementBases(basesEnd);
                suppBases = basesEnd + Nucleotides.reverseComplementBases(basesStart);
            }
            else
            {
                readBases = Nucleotides.reverseComplementBases(basesEnd) + basesStart;
                suppBases = Nucleotides.reverseComplementBases(basesStart) + basesEnd;
            }
        }

        List<SAMRecord> reads = Lists.newArrayList();

        String mateCigar = NO_CIGAR;
        int mateEnd = 0;
        String mateBases = null;

        if(mateStart > 0)
        {
            mateCigar = format("%dM", readBaseLength);
            mateEnd = mateStart + readBaseLength - 1;
            mateBases = refGenome.getBaseString(chrStart, mateStart, mateEnd);
        }

        SupplementaryReadData readSuppData = new SupplementaryReadData(
                chrEnd, suppStart, isSuppNegStrand ? SUPP_NEG_STRAND : SUPP_POS_STRAND, suppCigar, 60, 0);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readId, chrStart, readStart, readBases, readCigar, chrStart, mateStart, false,
                false, readSuppData, true, mateCigar);

        reads.add(read);

        SupplementaryReadData suppReadData = new SupplementaryReadData(
                chrStart, readStart, SUPP_POS_STRAND, readCigar, 60, 0);

        SAMRecord supp = SamRecordTestUtils.createSamRecord(
                readId, chrEnd, suppStart, suppBases, suppCigar, chrStart, mateStart, false,
                true, suppReadData, true, mateCigar);

        reads.add(supp);

        if(mateStart > 0)
        {
            SAMRecord mate = SamRecordTestUtils.createSamRecord(
                readId, chrStart, mateStart, mateBases, mateCigar, chrStart, readStart, true,
                false, null, false, readCigar);

            reads.add(mate);
        }

        return reads;
    }

    public static String formTestRefSequence(final int length)
    {
        // tries to avoid long repeats
        int currentIndex = 0;
        int maxSegmentLength = 40;

        StringBuilder sb = new StringBuilder();
        int currentLength = 0;

        while(currentLength < length)
        {
            String nextSegment = REF_BASES_400.substring(currentIndex, currentIndex + maxSegmentLength)
                    + MockRefGenome.generateRandomBases(10);

            if(nextSegment.length() + currentLength > length)
                nextSegment = nextSegment.substring(0, length - currentLength);

            sb.append(nextSegment);
            currentLength += nextSegment.length();

            ++currentIndex;

            if(currentIndex + maxSegmentLength > REF_BASES_400.length())
                currentIndex = 0;
        }

        return sb.toString();
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
        SAMRecord record = createSamRecord(
                readId, chromosome, readStart, "", "100M",
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

    public static void loadRefGenomeBases(final MockRefGenome refGenome, final String testFilename)
    {
        List<String> lines = new BufferedReader(new InputStreamReader(
                TestUtils.class.getResourceAsStream(testFilename))).lines().collect(Collectors.toList());

        for(String line : lines)
        {
            String[] values = line.split(CSV_DELIM, 2);
            String chr = values[0];

            if(chr.startsWith("#")) // comment lines
                continue;

            String bases = values[1];

            String existingBases = refGenome.RefGenomeMap.get(chr);

            if(existingBases != null)
                bases = existingBases + bases;

            refGenome.RefGenomeMap.put(chr, bases);
        }
    }
}