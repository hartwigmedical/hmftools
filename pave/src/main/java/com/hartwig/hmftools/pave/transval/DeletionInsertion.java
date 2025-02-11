package com.hartwig.hmftools.pave.transval;

import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;

class DeletionInsertion extends ProteinVariant
{
    @NotNull
    final private AminoAcidSequence Alt;

    DeletionInsertion(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            @NotNull final AminoAcidRange refRange,
            @NotNull final AminoAcidSequence alt)
    {
        super(gene, transcript, aminoAcidSequence, refRange.startPosition(), refRange.length());
        this.Alt = alt;
    }

    @VisibleForTesting
    public SplitCodonSequence referenceBases(RefGenomeInterface genome)
    {
        ChangeContext changeContext = getChangeContext(genome);
        return changeContext.basesForProteinChange(positionOfFirstAlteredCodon(), RefLength, Transcript.posStrand());
    }

    @Override
    public TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        SplitCodonSequence splitBases = referenceBases(refGenome);
        String referenceNucleotides = splitBases.completeSequence();
        if(!splitBases.couldBeDeletionInsertion())
        {
            return null; // todo
        }
        Set<String> destinationBases = candidateAlternativeNucleotideSequences(splitBases.retainedPrefix(), splitBases.retainedSuffix());
        String modifiedBases = Gene.forwardStrand()
                ? splitBases.segmentThatIsModified()
                : Nucleotides.reverseComplementBases(splitBases.segmentThatIsModified());
        int changeStart = splitBases.locationOfDeletedBases();
        if(Gene.reverseStrand())
        {
            destinationBases = destinationBases.stream().map(Nucleotides::reverseComplementBases).collect(Collectors.toSet());
            changeStart = changeStart - modifiedBases.length() + 1;
            referenceNucleotides = Nucleotides.reverseComplementBases(referenceNucleotides);
        }
        SortedSet<DeletionInsertionChange> changes = new TreeSet<>();
        destinationBases.forEach(target -> changes.add(new DeletionInsertionChange(modifiedBases, target)));

        ChangeLocation globalLocation = new ChangeLocation(this.Gene.Chromosome, changeStart);
        Set<TransvalHotspot> hotspots = changes
                .stream()
                .map(change -> change.toHotspot(globalLocation))
                .collect(Collectors.toSet());
        SortedSet<TransvalHotspot> hotspotsSorted = new TreeSet<>(hotspots);
        TransvalHotspot bestCandidate = hotspotsSorted.first();

        return new TransvalDeletionInsertion(
                Transcript,
                Gene.Chromosome,
                bestCandidate.mPosition,
                splitBases.spansTwoExons(),
                referenceNucleotides,
                bestCandidate.Ref,
                hotspots
        );
    }

    public String altAminoAcidSequence()
    {
        return Alt.sequence();
    }

    @VisibleForTesting
    public Set<String> candidateAlternativeNucleotideSequences(@NotNull String prefix, @NotNull String suffix)
    {
        return new NucleotidesCalculator(Alt, prefix, suffix).possibilities();
    }
}
