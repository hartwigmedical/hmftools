package com.hartwig.hmftools.lilac.misc;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_BASE_QUAL;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.Read;
import com.hartwig.hmftools.lilac.seq.HlaSequence;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

import htsjdk.samtools.SAMRecord;

public class LilacTestUtils
{
    public static void disableLogging()
    {
        Configurator.setRootLevel(Level.ERROR);
    }

    public static HlaSequenceLoci createSequenceLoci(final HlaSequence sequences)
    {
        return HlaSequenceLoci.create(sequences.Allele, sequences.getRawSequence(), sequences.getRawSequence());
    }

    public static final String TEST_READ_ID = "READ_ID_001";
    public static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(151);
    public static final String TEST_READ_CIGAR = "151M";

    public static Read createReadRecord(final String readId)
    {
        SAMRecord record = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 300,
                false, false, null);

        return Read.createRead(new BaseRegion(0, 1), record, true, true);
    }

    public static Fragment createFragment(final String id)
    {
        return new Fragment(
                createReadRecord(id), HlaGene.NONE, Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());
    }

    public static Fragment createFragment(final String id, final HlaGene gene, final String sequence, int locusStart, int locusEnd)
    {
        List<Integer> loci = formRange(locusStart, locusEnd);
        List<String> sequences = buildTargetSequences(sequence, loci);
        List<Byte> qualities = loci.stream().map(x -> DEFAULT_MIN_BASE_QUAL).collect(Collectors.toList());

        return new Fragment(createReadRecord(id), gene, Sets.newHashSet(gene), loci, qualities, sequences);
    }

    public static String buildTargetSequence(final String sequence, final List<Integer> indices)
    {
        if(indices.size() >= sequence.length())
            return "";

        StringBuilder targetSeq = new StringBuilder();
        indices.forEach(x -> targetSeq.append(sequence.charAt(x)));
        return targetSeq.toString();
    }

    public static List<Integer> buildLoci(final String sequence)
    {
        List<Integer> indices = Lists.newArrayList();

        for(int i = 0; i < sequence.length(); ++i)
        {
            indices.add(i);
        }

        return indices;
    }


    public static List<String> buildTargetSequences(final String sequence, final List<Integer> indices)
    {
        List<String> sequences = Lists.newArrayList();
        if(indices.size() > sequence.length())
            return sequences;

        indices.forEach(x -> sequences.add(String.valueOf(sequence.charAt(x))));
        return sequences;
    }

    public static List<Integer> buildTargetLoci(final String sequence, final String refSequence)
    {
        List<Integer> indices = Lists.newArrayList();

        if(sequence.length() != refSequence.length())
            return indices;

        for(int i = 0; i < sequence.length(); ++i)
        {
            if(sequence.charAt(i) != refSequence.charAt(i))
            {
                indices.add(i);
            }
        }

        return indices;
    }
}
