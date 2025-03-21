package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.util.Checks;

import org.apache.commons.lang3.tuple.Pair;

public class SubstitutionVariant extends DnaVariant
{
    public final String Ref;
    public final String Alt;

    public SubstitutionVariant(GeneData gene, TranscriptData transcript, HgvsAddress start, HgvsAddress finish, String ref, String alt)
    {
        super(gene, transcript, start, finish);
        Preconditions.checkArgument(Checks.isNucleotideSequence(ref));
        Preconditions.checkArgument(Checks.isNucleotideSequence(alt));
        Ref = ref;
        Alt = alt;
    }

    public SubstitutionVariant(GeneData gene, TranscriptData transcript, HgvsAddress address, String ref, String alt)
    {
        this(gene, transcript, address, address, ref, alt);
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        Pair<Integer, Integer> startStop = getAbsoluteLocationsOfChange();
        String refFromGenome = genome.getBaseString(chromosome(), startStop.getLeft(), startStop.getRight());
        String alt = forwardStrand() ? Alt : Nucleotides.reverseComplementBases(Alt);
        return new BaseSequenceChange(refFromGenome, alt, chromosome(), startStop.getLeft());
    }
}
