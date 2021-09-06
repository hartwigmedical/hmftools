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

    public static boolean isFullyGTelomeric(final String readBases)
    {
        if (StringUtils.countMatches(readBases, "TTAGGG") >= MIN_CANONICAL_COUNT)
        {
            Matcher m = sGTeloPatternMatcher.get();
            m.reset(readBases);
            return m.find();
        }
        return false;
    }

    public static boolean isFullyCTelomeric(final String readBases)
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
        boolean isGTelomeric = isFullyGTelomeric(readPairG);
        boolean isCTelomeric = isFullyCTelomeric(readPairC);

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
        if(!config.RefGenomeFile.isEmpty())
        {
            factory = factory.referenceSequence(new File(config.RefGenomeFile));
        }
        return factory.open(new File(config.BamFile));
    }
}
