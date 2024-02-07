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
import static com.hartwig.hmftools.esvee.read.ReadUtils.addByteArray;
import static com.hartwig.hmftools.esvee.read.ReadUtils.reverseBytes;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.SAMRecord;

public class TestUtils
{
    public static final String REF_BASES = generateRandomBases(100);
    public static final String TEST_READ_ID = "READ_01";
    public static final String TEST_CIGAR  = "100M";

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
            final String chrEnd, int junctionPosEnd, byte junctionOrientEnd, int matePosStart)
    {
        int readBaseLength = anchorLength * 2;
        int readStart, readEnd, suppStart, suppEnd, mateStart;
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
                readBases = basesStart + Nucleotides.reverseStrandBases(basesEnd);
                suppBases = basesEnd + Nucleotides.reverseStrandBases(basesStart);
            }
            else
            {
                readBases = Nucleotides.reverseStrandBases(basesEnd) + basesStart;
                suppBases = Nucleotides.reverseStrandBases(basesStart) + basesEnd;
            }
        }

        List<SAMRecord> reads = Lists.newArrayList();

        String mateCigar = NO_CIGAR;
        int matePosEnd = 0;
        byte[] mateBases = null;

        if(matePosStart > 0)
        {
            mateCigar = format("%dM", anchorLength * 2);
            matePosEnd = matePosStart + anchorLength * 2 - 1;
            mateBases = refGenome.getBases(chrStart, matePosStart, matePosEnd);
        }

        SupplementaryReadData readSuppData = new SupplementaryReadData(
                chrEnd, suppStart, isSuppNegStrand ? SUPP_NEG_STRAND : SUPP_POS_STRAND, suppCigar, 60, 0);

        SupplementaryReadData suppReadData = new SupplementaryReadData(
                chrStart, readStart, SUPP_POS_STRAND, readCigar, 60, 0);

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readId, chrStart, readStart, readBases, readCigar, chrStart, matePosStart, false,
                false, readSuppData, true, mateCigar);

        reads.add(read);

        /*
        final String readId, final String chrStr, int readStart, final String readBases, final String cigar, final String mateChr,
            int mateStart, boolean isReversed, boolean isSupplementary, final SupplementaryReadData suppAlignment,
            boolean mateReversed, final String mateCigar
        */


        return reads;
    }
}