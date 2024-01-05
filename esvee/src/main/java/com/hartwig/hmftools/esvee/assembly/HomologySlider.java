package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.sequence.AlignedAssembly;
import com.hartwig.hmftools.esvee.sequence.Alignment;

/** This slides everything left for deduplication purposes, we will slide it back to the center in a later step */
public class HomologySlider
{
    private final RefGenomeInterface mReference;

    // Input Alignment / Sequence
    // Left                 Right
    // 01234567890123456789 01234567890123456789
    // AAAAAAAAAAAAAAAAAAAA AAAATTTTTTTTTTTTTTTT
    // Sequence Alignment
    // AAAAAAAAAAAAAAAA     AAAATTTTTT
    // Desired Alignment
    // Left                 Right
    // 01234567890123456789 01234567890123456789
    // AAAAAAAAAAAAAAAAAAAA AAAATTTTTTTTTTTTTTTT
    // Sequence Alignment
    // AAAAAAAAAAAAAAAAAAAA     TTTTTT

    public HomologySlider(final RefGenomeInterface reference)
    {
        mReference = reference;
    }

    @SuppressWarnings("ConstantValue")
    public List<Alignment> slideHomology(final List<Alignment> source)
    {
        // For each part of the alignment, keep expanding the LHS if possible
        final List<Alignment> alignments = new ArrayList<>();
        alignments.add(source.get(0));
        for(int i = 1; i < source.size(); i++)
        {
            final Alignment left = alignments.get(alignments.size() - 1);
            final Alignment right = source.get(i);
            if(left.Inverted || right.Inverted)
                return source; // FIXME: Not handled yet
            if(left.isUnmapped() || right.isUnmapped())
            {
                alignments.add(right);
                continue;
            }

            final int leftEnd = left.ReferenceStartPosition + left.Length;
            final int leftMax = Math.min(leftEnd + 100, mReference.getChromosomeLength(left.Chromosome));
            final byte[] afterLeftBases = mReference.getBases(left.Chromosome, leftEnd, leftMax);
            final byte[] rightBases =
                    mReference.getBases(right.Chromosome, right.ReferenceStartPosition, right.ReferenceStartPosition + right.Length - 1);

            int homology = 0;
            for(int j = 0; j < Math.min(afterLeftBases.length, rightBases.length); j++)
                if(afterLeftBases[j] == rightBases[j])
                    homology++;
                else
                    break;

            if(homology > 0)
            {
                // Expand left to cover right
                final Alignment oldLeft = alignments.get(alignments.size() - 1);
                final Alignment newLeft = new Alignment(oldLeft.Chromosome, oldLeft.ReferenceStartPosition,
                        oldLeft.SequenceStartPosition, oldLeft.Length + homology, oldLeft.Inverted, oldLeft.Quality);
                alignments.set(alignments.size() - 1, newLeft);
                final int newLength = right.Length - homology;
                if(newLength > 0)
                {
                    final Alignment newRight = new Alignment(right.Chromosome, right.ReferenceStartPosition + homology,
                            right.SequenceStartPosition + homology, newLength, right.Inverted, right.Quality);
                    alignments.add(newRight);
                }
            }
            else
                alignments.add(right);
        }

        return alignments;
    }

    public AlignedAssembly slideHomology(final AlignedAssembly assembly)
    {
        final AlignedAssembly aligned = new AlignedAssembly(assembly.Source, slideHomology(assembly.getAlignmentBlocks()));
        return aligned;
    }
}
