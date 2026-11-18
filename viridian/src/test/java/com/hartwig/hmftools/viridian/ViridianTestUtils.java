package com.hartwig.hmftools.viridian;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusReference;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;

public final class ViridianTestUtils
{
    public static ViralReadAlignment readAlignment(
            String readName, VirusGenome genome, int alignmentStart, String cigar, int alignerScore, int divergence)
    {
        return new ViralReadAlignment(
                ReadId.parse(readName), genome, alignmentStart, TextCigarCodec.decode(cigar), alignerScore, divergence,
                new int[0]);
    }

    public static VirusReference virusReference(VirusGenome... genomes)
    {
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(Arrays.stream(genomes)
                .map(genome -> new SAMSequenceRecord(genome.contigName(), genome.length()))
                .toList());
        return new VirusReference(List.of(genomes), dictionary, Map.of());
    }

    public static int samFlags(SAMFlag... flags)
    {
        return Arrays.stream(flags).mapToInt(SAMFlag::intValue).reduce(0, (a, b) -> a | b);
    }

    // What BWA-MEM returns for a query that aligns nowhere.
    public static BwaMemAlignment noHit()
    {
        return new BwaMemAlignment(samFlags(SAMFlag.READ_UNMAPPED), -1, -1, -1, 0, 0, 0, 0, 0, 0, "*", null, null, -1, -1, 0);
    }
}
