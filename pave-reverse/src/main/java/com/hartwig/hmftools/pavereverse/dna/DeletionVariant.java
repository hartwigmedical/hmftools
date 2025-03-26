package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.apache.commons.lang3.tuple.Pair;

public class DeletionVariant extends DnaVariant
{
    private final String mDeletedBases;

    public DeletionVariant(GeneData gene, TranscriptData transcript, HgvsAddress address, String deletedBases)
    {
        this(gene, transcript, address, address, deletedBases);
    }

    public DeletionVariant(GeneData gene, TranscriptData transcript, HgvsAddress start, HgvsAddress end, String deletedBases)
    {
        super(gene, transcript, start, end);
        mDeletedBases = deletedBases;
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        Pair<Integer, Integer> startStop = getAbsoluteLocationsOfChange();
        LeftMostEquivalentDeletionFinder finder = new LeftMostEquivalentDeletionFinder(genome, chromosome(), startStop.getLeft(), startStop.getRight());
        int leftMostPosition = finder.findLeftMostEquivalentPosition();
        int neighbourPosition = leftMostPosition - 1;
        int stop = leftMostPosition + startStop.getRight() - startStop.getLeft();
        String ref = genome.getBaseString(chromosome(), neighbourPosition, stop);
        return new BaseSequenceChange(ref, ref.substring(0, 1), chromosome(), neighbourPosition);
    }
}
