package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConstants.*;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReader;
import org.apache.commons.lang3.StringUtils;

public class TeloUtils
{
    static final List<String> CANONICAL_TELOMERE_SEQUENCES = com.google.common.collect.Lists.newArrayList(
            String.join("", Collections.nCopies(DEFAULT_MIN_TELE_SEQ_COUNT, CANONICAL_TELOMERE_SEQ)),
            String.join("", Collections.nCopies(DEFAULT_MIN_TELE_SEQ_COUNT, CANONICAL_TELOMERE_SEQ_REV)));

    public static boolean hasTelomericContent(final String readBases)
    {
        return hasTelomericContent(readBases, CANONICAL_TELOMERE_SEQUENCES);
    }

    public static boolean hasTelomericContent(final String readBases, final List<String> sequences)
    {
        for(String teloSeq : sequences)
        {
            int matchIndex = readBases.indexOf(teloSeq);

            if (matchIndex != -1)
            {
                return true;
            }
        }

        return false;
    }

    // todo: try SequenceUtil.reverseComplement
    public static String reverseComplementSequence(String seq)
    {
        StringBuilder builder = new StringBuilder();

        for(int i = seq.length() - 1; i >= 0; i--)
        {
            char b = seq.charAt(i);
            if (b == 'A')
            {
                builder.append('T');
            }
            else if (b == 'T')
            {
                builder.append('A');
            }
            else if (b == 'G')
            {
                builder.append('C');
            }
            else if (b == 'C')
            {
                builder.append('G');
            }
        }

        return builder.toString();
    }

    private static final int MIN_CANONICAL_COUNT = 4;
    private static final int MIN_CONSECUTIVE_HEXAMERS = 6;

    // we match the start of the sequence, must be at least 6 repeats of telomere. The reason we only match the start of the sequence
    // is to account for poly G tail that many reads have
    private static Pattern sGTeloPattern = Pattern.compile("^[ACGT]{0,2}G{0,3}" + StringUtils.repeat("(?!GGGGGG)([ACGT][ACGT][ACGT]GGG)", MIN_CONSECUTIVE_HEXAMERS));
    private static Pattern sCTeloPattern = Pattern.compile("^[ACGT]{0,2}C{0,3}" + StringUtils.repeat("(?!CCCCCC)([ACGT][ACGT][ACGT]CCC)", MIN_CONSECUTIVE_HEXAMERS));

    // make matcher thread local since they are not thread safe
    private static ThreadLocal<Matcher> sGTeloPatternMatcher = ThreadLocal.withInitial(() -> sGTeloPattern.matcher(""));
    private static ThreadLocal<Matcher> sCTeloPatternMatcher = ThreadLocal.withInitial(() -> sCTeloPattern.matcher(""));

    public static boolean isLikelyGTelomeric(final String readBases)
    {
        if (StringUtils.countMatches(readBases, "TTAGGG") >= MIN_CANONICAL_COUNT)
        {
            Matcher m = sGTeloPatternMatcher.get();
            m.reset(readBases);
            return m.find();
        }
        return false;
    }

    public static boolean isLikelyCTelomeric(final String readBases)
    {
        if (StringUtils.countMatches(readBases, "CCCTAA") >= MIN_CANONICAL_COUNT)
        {
            Matcher m = sCTeloPatternMatcher.get();
            m.reset(readBases);
            return m.find();
        }
        return false;
    }

    public static ReadGroup.FragmentType classifyFragment(final ReadGroup readGroup)
    {
        if (readGroup.Reads.size() < 2)
        {
            return ReadGroup.FragmentType.UNKNOWN;
        }

        // we want to check
        SAMRecord read1 = readGroup.Reads.get(0);
        SAMRecord read2 = readGroup.Reads.get(1);

        String seq1 = read1.getReadString();
        String seq2 = read2.getReadString();

        // qual_str1 = read1['BaseQualities']
        // qual_str2 = read2['BaseQualities']

        if (read1.getReadNegativeStrandFlag())
        {
            seq1 = reverseComplementSequence(seq1);
        }

        if (read2.getReadNegativeStrandFlag())
        {
            seq2 = reverseComplementSequence(seq2);
        }

        // decide which one is G seq and which one is C seq, by counting GGG
        if (StringUtils.countMatches(seq1, "TTAGGG")  > StringUtils.countMatches(seq2, "TTAGGG"))
        {
            return classifyFragment(seq1, seq2);
        }
        else
        {
            return classifyFragment(seq2, seq1);
        }
    }

    public static ReadGroup.FragmentType classifyFragment(String readPairG, String readPairC)
    {
        boolean isGTelomeric = isLikelyGTelomeric(readPairG);
        boolean isCTelomeric = isLikelyCTelomeric(readPairC);

        // first we want to find the TTAGGG motif, then fill it up backwards
        if (isGTelomeric && isCTelomeric)
        {
            return ReadGroup.FragmentType.F1;
        }

        if (isGTelomeric)
        {
            // only g is telomeric
            return ReadGroup.FragmentType.F4;
        }

        if (isCTelomeric)
        {
            // only C term is telomeric, could be F2a or F2b
            return ReadGroup.FragmentType.F2;
        }

        // not telomeric
        return ReadGroup.FragmentType.NOT_TELOMERE;
    }

    // TTAGGG but we also match for anything that start with [ACGT]{3}GGG with total <= 4 Gs
    static boolean isGTeloHexamer(String sequence, int startOffset)
    {
        if ((sequence.length() - startOffset) < 6)
        {
            return false;
        }

        for (String gHexamer : TELOMERE_HEXAMERS)
        {
            if (sequence.regionMatches(startOffset, gHexamer, 0, gHexamer.length()))
                return true;
        }

        return false;
    }

    // CCCTAA but we also match for anything that start with CCC[ACGT]{3} with total <= 4 Cs
    static boolean isCTeloHexamer(String sequence, int startOffset)
    {
        if ((sequence.length() - startOffset) < 6)
        {
            return false;
        }

        for (String cHexamer : TELOMERE_HEXAMERS_REV)
        {
            if (sequence.regionMatches(startOffset, cHexamer, 0, cHexamer.length()))
                return true;
        }

        return false;
    }

    public static boolean isStrictlyGTelomeric(final String seq)
    {
        return isStrictlyTelomericHelper(seq, CANONICAL_TELOMERE_SEQ, TeloUtils::isGTeloHexamer);
    }

    public static boolean isStrictlyCTelomeric(final String seq)
    {
        return isStrictlyTelomericHelper(seq, CANONICAL_TELOMERE_SEQ_REV, TeloUtils::isCTeloHexamer);
    }

    public static boolean isStrictlyTelomericHelper(final String seq, final String canonicalHexamer,
            java.util.function.BiFunction<String, Integer, Boolean> hexamerCheckFunc)
    {
        // here we apply a very stringent rule to go through the bases and work out if something is
        // telomeric by shifting across by 6 and see if we can find a way that is totally telomeric
        for (int i = 0; i < 6; ++i)
        {
            // we use a modified seq to fix up the ends to make it easier
            String modifiedSeq = seq;
            if (i != 0)
            {
                // we add something at the start to make it a full hexamer
                modifiedSeq = canonicalHexamer.substring(0, 6 - i) + seq;
            }

            // also fix the end
            int residualLength = modifiedSeq.length() % 6;
            if (residualLength != 0)
            {
                modifiedSeq += canonicalHexamer.substring(residualLength);
            }

            boolean isAllTelo = true;
            for (int j = 0; j < modifiedSeq.length(); j += 6)
            {
                if (!hexamerCheckFunc.apply(modifiedSeq, j))
                {
                    isAllTelo = false;
                    break;
                }
            }

            if (isAllTelo)
                return true;
        }
        return false;
    }

    // determine if the read is poly G
    // note: We must supply this function with the read as is from
    // the read
    public static boolean isPolyGC(String readSeq)
    {
        // we use threshold of 90%
        double numGCs = Math.max(StringUtils.countMatches(readSeq, 'G'), StringUtils.countMatches(readSeq, 'C'));
        double gcFrac = numGCs / readSeq.length();
        return gcFrac >= POLY_G_THRESHOLD;
    }

    public static List<ChrBaseRegion> createPartitions(final TeloConfig config)
    {
        SamReader samReader = TeloUtils.openSamReader(config);

        List<SAMSequenceRecord> samSequences = samReader.getFileHeader().getSequenceDictionary().getSequences();

        List<ChrBaseRegion> partitions = Lists.newArrayList();

        int partitionSize = DEFAULT_PARTITION_SIZE;

        for(SAMSequenceRecord seq : samSequences)
        {
            String chrStr = seq.getSequenceName();

            if(!config.SpecificChromosomes.isEmpty() && !config.SpecificChromosomes.contains(chrStr))
                continue;

            int chromosomeLength = seq.getSequenceLength();

            int startPos = 0;
            while(startPos < chromosomeLength)
            {
                int endPos = startPos + partitionSize - 1;

                if(endPos + partitionSize * 0.2 > chromosomeLength)
                    endPos = chromosomeLength;

                partitions.add(new ChrBaseRegion(chrStr, startPos, endPos));

                startPos = endPos + 1;
            }
        }

        return partitions;
    }

    public static SamReader openSamReader(final TeloConfig config)
    {
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        if(config.RefGenomeFile != null && !config.RefGenomeFile.isEmpty())
        {
            factory = factory.referenceSequence(new File(config.RefGenomeFile));
        }
        return factory.open(new File(config.BamFile));
    }
}
