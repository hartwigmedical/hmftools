package com.hartwig.hmftools.pave.transval;

import static com.hartwig.hmftools.pave.transval.Checks.isValidProtein;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class DeletionInsertion extends ProteinVariant
{
    final private int DeletionLength;
    @NotNull
    final private String Alt;

    public DeletionInsertion(
            @NotNull final GeneData gene,
            @NotNull final TranscriptData transcript,
            @NotNull final TranscriptAminoAcids aminoAcidSequence,
            final int position,
            final int deletionLength,
            @NotNull final String variant)
    {
        super(gene, transcript, aminoAcidSequence, position);
        Preconditions.checkArgument(isValidProtein(variant));
        this.Alt = variant;
        this.DeletionLength = deletionLength;
    }

    @VisibleForTesting
    public SplitSequence referenceBases(RefGenomeInterface genome)
    {
        int codonPosition = 3 * (positionOfFirstAlteredCodon() - 1);
        List<Integer> regionLengths = codingRegionLengths();
        int lengthIncludingCurrent = 0;
        for(int i = 0; i < regionLengths.size(); i++)
        {
            int lengthUpToCurrent = lengthIncludingCurrent;
            lengthIncludingCurrent += regionLengths.get(i);
            ChrBaseRegion exon = CodingRegions.get(i);
            if(lengthIncludingCurrent > codonPosition)
            {
                ChrBaseRegion nextExon = (i < CodingRegions.size() - 1) ? CodingRegions.get(i + 1) : null;
                if(Transcript.negStrand())
                {
                    //                        int positionOfCodonInCurrentExon = exon.end() - (codonPosition - lengthUpToCurrent);
                    //                        return new CodonRegions(positionOfCodonInCurrentExon, exon, nextExon, false);
                    return null;
                }
                else
                {
                    int positionOfStartInCurrent = exon.start() + (codonPosition - lengthUpToCurrent);
                    int positionOfEnd = positionOfStartInCurrent + DeletionLength * 3 - 1;
                    if(positionOfEnd <= exon.end())
                    {
                        String partInCurrent = genome.getBaseString(Gene.Chromosome, positionOfStartInCurrent, positionOfEnd);
                        return new SplitSequence(partInCurrent, null, positionOfStartInCurrent);
                    }
                    else
                    {
                        if(nextExon == null)
                        {
                            throw new IllegalStateException("Variant does not fit given exons.");
                        }
                        String partInCurrentExon = genome.getBaseString(Gene.Chromosome, positionOfStartInCurrent, exon.end());
                        int lengthInNext = positionOfEnd - exon.end();
                        // If the right hand part is of length 3 or more, then the
                        // changes are on that side only, and the left hand side is
                        // a fixed prefix.
                        int positionOfChange = lengthInNext > 2 ? nextExon.start() : positionOfStartInCurrent;
                        String partInNext = genome.getBaseString(Gene.Chromosome, nextExon.start(), nextExon.start() + lengthInNext - 1);
                        return new SplitSequence(partInCurrentExon, partInNext, positionOfChange);
                    }
                }
            }
        }
        //            throw new IllegalArgumentException("No exon found for codon " + codonPosition);

        return null;
    }

    @Override
    public TransvalVariant calculateVariant(RefGenomeInterface refGenome)
    {
        SplitSequence splitBases = referenceBases(refGenome);
        if(!splitBases.couldBeDeletionInsertion())
        {
            return null; // todo
        }
        Set<String> destinationBases = candidateAlternativeNucleotideSequences(splitBases.retainedPrefix(), splitBases.retainedSuffix());
        String modifiedBases = splitBases.segmentThatIsModified();
        SortedSet<DeletionInsertionChange> changes = new TreeSet<>();
        destinationBases.forEach(target -> changes.add(new DeletionInsertionChange(modifiedBases, target)));

        ChangeLocation globalLocation = new ChangeLocation(this.Gene.Chromosome, splitBases.locationOfDeletedBases());
        Set<TransvalHotspot> hotspots = changes
                .stream()
                .map(change -> change.toHotspot(globalLocation))
                .collect(Collectors.toSet());
        DeletionInsertionChange bestCandidate = changes.first();
        int positionOfDeletionStart = splitBases.locationOfDeletedBases() + bestCandidate.positionOfDeletion();

        return new TransvalInsertionDeletion(
                Transcript,
                Gene.Chromosome,
                positionOfDeletionStart,
                splitBases.spansTwoExons(),
                splitBases.completeSequence(),
                bestCandidate.deleted(),
                hotspots
        );
    }

    public String altAminoAcidSequence()
    {
        return Alt;
    }

    @Override
    int changedReferenceSequenceLength()
    {
        return DeletionLength;
    }

    @VisibleForTesting
    public Set<String> candidateAlternativeNucleotideSequences(@NotNull String prefix, @NotNull String suffix)
    {
        return new NucleotidesCalculator(Alt, prefix, suffix).possibilities();
    }
}
