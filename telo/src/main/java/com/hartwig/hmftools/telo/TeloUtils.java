package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConstants.*;

import java.io.File;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReader;

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

    public static List<BaseRegion> createPartitions(final TeloConfig config)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.RefGenomeFile)).open(new File(config.SampleBamFile));

        List<SAMSequenceRecord> samSequences = samReader.getFileHeader().getSequenceDictionary().getSequences();

        List<BaseRegion> partitions = Lists.newArrayList();

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

                partitions.add(new BaseRegion(chrStr, startPos, endPos));

                startPos = endPos + 1;
            }
        }

        return partitions;
    }
}
