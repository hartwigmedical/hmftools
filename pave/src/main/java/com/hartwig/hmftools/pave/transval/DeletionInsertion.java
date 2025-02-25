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
        return changeContext.basesForProteinChange(positionOfFirstAlteredCodon(), mRefLength, mTranscript.posStrand());
    }

    @Override
    public TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        SplitCodonSequence splitBases = referenceBases(refGenome);
        String referenceNucleotidesForwardStrand = splitBases.completeSequence();
        if(mGene.reverseStrand())
        {
            referenceNucleotidesForwardStrand = Nucleotides.reverseComplementBases(referenceNucleotidesForwardStrand);
        }
        if(!splitBases.couldBeDeletionInsertion())
        {
            return null; // todo
        }
        Set<String> destinationBasesForwardStrand = candidateAlternativeNucleotideSequences(splitBases.retainedPrefix(), splitBases.retainedSuffix());
        if(mGene.reverseStrand())
        {
            destinationBasesForwardStrand = destinationBasesForwardStrand.stream().map(Nucleotides::reverseComplementBases).collect(Collectors.toSet());
            //            changeStart = changeStart - modifiedForwardStrandBases.length() + 1;
        }
        int changeStart = splitBases.locationOfDeletedBases();
        String modifiedForwardStrandBases = mGene.forwardStrand()
                ? splitBases.segmentThatIsModified()
                : Nucleotides.reverseComplementBases(splitBases.segmentThatIsModified());
        SortedSet<DeletionInsertionChange> changes = new TreeSet<>();
        destinationBasesForwardStrand.forEach(target -> changes.add(new DeletionInsertionChange(modifiedForwardStrandBases, target)));

        ChangeLocation globalLocation = new ChangeLocation(this.mGene.Chromosome, changeStart);
        Set<TransvalHotspot> hotspots = changes
                .stream()
                .map(change -> change.toHotspot(globalLocation))
                .collect(Collectors.toSet());
        SortedSet<TransvalHotspot> hotspotsSorted = new TreeSet<>(hotspots);
        TransvalHotspot bestCandidate = hotspotsSorted.first();

        return new TransvalDeletionInsertion(
                mTranscript,
                mGene.Chromosome,
                bestCandidate.mPosition,
                splitBases.spansTwoExons(),
                referenceNucleotidesForwardStrand,
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
        return new NucleotidesCalculator(Alt, prefix, suffix).allPossibleBaseSequences();
    }
}
