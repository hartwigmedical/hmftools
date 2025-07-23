package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.apache.commons.lang3.tuple.Pair;

public class InsertionVariant extends DnaVariant
{
    private final String mInsertedBases;

    public InsertionVariant(GeneData gene, TranscriptData transcript, HgvsAddress start, HgvsAddress end, String insertedBases)
    {
        super(gene, transcript, start, end);
        mInsertedBases = insertedBases;
    }

    public String insertedBases()
    {
        return mInsertedBases;
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        Pair<Integer, Integer> startStop = getAbsoluteLocationsOfChange();
        String insertedBases = forwardStrand() ? mInsertedBases : Nucleotides.reverseComplementBases(mInsertedBases);
        LeftMostEquivalentInsertionFinder finder =
                new LeftMostEquivalentInsertionFinder(genome, chromosome(), startStop.getLeft(), insertedBases);
        int canonicalInsertionPoint = finder.findLeftMostEquivalentPosition();
        String ref = genome.getBaseString(chromosome(), canonicalInsertionPoint, canonicalInsertionPoint);
        return new BaseSequenceChange(ref, ref + insertedBases, chromosome(), canonicalInsertionPoint);
    }
}
