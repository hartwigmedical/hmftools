package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.apache.commons.lang3.tuple.Pair;

public class DuplicationVariant extends DnaVariant
{
    private final String mDuplicatedBases;

    public DuplicationVariant(GeneData gene, TranscriptData transcript, HgvsAddress start, HgvsAddress end, String duplicatedBases)
    {
        super(gene, transcript, start, end);
        mDuplicatedBases = duplicatedBases;
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        Pair<Integer, Integer> startStop = getAbsoluteLocationsOfChange();
        int canonicalStart = new LeftMostEquivalentDuplicationFinder(genome, chromosome(), startStop.getLeft(), startStop.getRight()).findLeftMostEquivalentPosition();
        int neighbourPosition = canonicalStart - 1;
        int end = neighbourPosition + startStop.getRight() - startStop.getLeft() + 1;
        String alt = genome.getBaseString(chromosome(), neighbourPosition, end);

        return new BaseSequenceChange(alt.substring(0, 1), alt, chromosome(), neighbourPosition);
    }
}
