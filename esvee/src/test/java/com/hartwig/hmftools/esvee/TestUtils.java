package com.hartwig.hmftools.esvee;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.samtools.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.cloneSamRecord;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.SAMRecord;

public class TestUtils
{
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_READ_ID_2 = "READ_02";

    public static final String TEST_CIGAR_100  = "100M";
    public static final String TEST_CIGAR_50  = "50M";
    public static final String TEST_CIGAR_30  = "30M";
    public static final String TEST_CIGAR_20  = "20M";

    public static final String REF_BASES_RANDOM_100 = generateRandomBases(100);

    public static String REF_BASES_200 =
            "AAACCCGGGTTTACGTAACCGGTTACGTAAAAACCCCCGGGGGTTTTTACGTAACCGGTTACGTAAACCCGGGTTTAAACGTTTTTGGGGCCCCAAAAAC"
    //       0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    // index 0         10        20        30        40        50        60        70        80        90
          + "GTGGCCTTAAACGTCCCAAAATTTTGGGGACGTGGGGGCCCCCAAAATTTTACGTCCGGTTAAACGTTTTCCCGGGAAAACGTAACCGGTTGGCCAAGCT";

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

    public static JunctionAssembly createAssembly(
            final String chromosome, final int junctionPosition, final byte junctionOrientation,
            final String assemblyBases, final int junctionIndex)
    {
        Junction junction = new Junction(chromosome, junctionPosition, junctionOrientation);

        int baseLength = assemblyBases.length();
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(baseLength);

        return new JunctionAssembly(junction, assemblyBases.getBytes(), baseQuals, junctionIndex);
    }

    public static Read createSamRecord(final String readId, int readStart, final String readBases, final String cigar)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, cigar, CHR_1, readStart + 1000,
                false, false, null);
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
            final String chrStart, int junctionPosStart, byte junctionOrientStart,
            final String chrEnd, int junctionPosEnd, byte junctionOrientEnd, int mateStart)
    {
        // creates a junction read, its supplementary and a local mate if the coords are supplied
        int readBaseLength = anchorLength * 2;
        int readStart, readEnd, suppStart, suppEnd;
        String readCigar, suppCigar;
        String basesStart, basesEnd;

        if(junctionOrientStart == POS_ORIENT)
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

        if(junctionOrientEnd == POS_ORIENT)
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
            if(junctionOrientStart == POS_ORIENT)
                readBases = basesStart + basesEnd;
            else
                readBases = basesEnd + basesStart;

            suppBases = readBases;
        }
        else
        {
            isSuppNegStrand = false;

            // keep the first read's bases in the 5' to 3' direction
            if(junctionOrientStart == POS_ORIENT)
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
}