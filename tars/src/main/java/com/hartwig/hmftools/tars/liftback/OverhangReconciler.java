package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.overhang.OverhangGate;

import htsjdk.samtools.SAMRecord;

// Applies the overhang gate to a record's lifted alignments: collapses each candidate, then finalizes the chosen
// primary and its tx-match XA alts. No-op when no gate is loaded.
// TODO revisit this later to simplify
public class OverhangReconciler
{
    private final OverhangGate mOverhangGate;

    public OverhangReconciler(final OverhangGate overhangGate)
    {
        mOverhangGate = overhangGate;
    }

    // Collapses each candidate's terminal overhangs against the genome; multiple candidates are also genomic-scored for
    // the Step 2 pick. Collapse only, no extension; a collapsed XA alt is dropped from the tag and locus count, never self.
    public void reconcileAlignmentsToGenome(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(mOverhangGate == null || !mOverhangGate.enabled() || alignments.isEmpty())
        {
            return;
        }

        byte[] forwardBases = record.getReadBases();
        if(forwardBases == null || forwardBases.length == 0)
        {
            return;
        }

        LiftedAlignment self = alignments.get(0);
        boolean recordForward = !record.getReadNegativeStrandFlag();
        byte[] reverseBases = null;

        // A single candidate has nothing to be ranked against, so it is not scored - and a mergeable supplementary
        // (same chromosome, within the max intron span) marks a genuine split read whose pre-merge halves must not be
        // scored either. A cross- or far-chromosome supplementary is scored as normal.
        boolean scoreCandidates = alignments.size() > 1 && !hasMergeableSupplementary(record);

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alignment = alignments.get(i);
            if(alignment.LiftedCigar == null)
            {
                continue;
            }

            // read bases in this candidate's orientation (reverse-complemented for an opposite-strand alt).
            byte[] bases;
            if(alignment.ForwardStrand == recordForward)
            {
                bases = forwardBases;
            }
            else
            {
                if(reverseBases == null)
                {
                    reverseBases = reverseComplement(forwardBases);
                }
                bases = reverseBases;
            }

            OverhangGate.Result collapsed = mOverhangGate.collapseJunctions(
                    alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, bases);

            // A collapsed XA alt is a fabricated placement: drop it from XA and the locus count. Self is never dropped.
            if(alignment != self && collapsed.dropped())
            {
                alignment.Dropped = true;
                continue;
            }

            if(collapsed.pos() != alignment.LiftedPos || !collapsed.cigar().equals(alignment.LiftedCigar))
            {
                alignment = alignment.withLiftedCigar(collapsed.pos(), collapsed.cigar());
                alignments.set(i, alignment);
            }

            // recompute the genomic score so candidates rank on one comparable scale (bwa's AlignmentScore is not
            // comparable across tx and ref)
            if(scoreCandidates)
            {
                alignment.GenomicScore = mOverhangGate.genomicScore(
                        alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, bases);
            }
        }
    }

    // True when the record has a supplementary on the same reference as the primary within the max intron span - a
    // genuine splice partner. Reads the supplementary's own SA-declared position.
    private static boolean hasMergeableSupplementary(final SAMRecord record)
    {
        List<SupplementaryReadData> supplementaries = SupplementaryReadData.extractAlignments(record);
        if(supplementaries == null)
        {
            return false;
        }

        String primaryChromosome = record.getReferenceName();
        int primaryPosition = record.getAlignmentStart();
        for(SupplementaryReadData supplementary : supplementaries)
        {
            if(supplementary.Chromosome.equals(primaryChromosome)
                    && Math.abs(supplementary.Position - primaryPosition) <= TarsConstants.MAX_IMPLIED_INTRON_LENGTH)
            {
                return true;
            }
        }
        return false;
    }

    // Finalizes the chosen primary: collapse terminal overhangs, then extend over a surviving terminal softclip when the
    // primary is a tx-match. A clean primary is unchanged.
    public void reconcileChosenPrimary(
            final LiftedRecord[] resolved, final SAMRecord primary, final int primaryIdx,
            final boolean[] droppedBySupplementary)
    {
        if(mOverhangGate == null || !mOverhangGate.enabled() || primary == null || primary.getReadUnmappedFlag())
        {
            return;
        }

        if(primaryIdx < 0 || (droppedBySupplementary != null && droppedBySupplementary[primaryIdx]))
        {
            return;
        }

        LiftedRecord primaryRes = resolved[primaryIdx];
        if(primaryRes == null || !primaryRes.hasPlacement())
        {
            return;
        }

        // reverse-complement the read bases when the resolved strand differs from the record's own, so the gate walks
        // them in the placement's orientation.
        byte[] readBases = primary.getReadBases();
        if(readBases != null && primaryRes.negativeStrand() != primary.getReadNegativeStrandFlag())
        {
            readBases = reverseComplement(readBases);
        }

        OverhangGate.Result collapsed = mOverhangGate.collapseJunctions(
                primaryRes.finalChromosome(), primaryRes.finalPos(), primaryRes.finalCigar(), readBases);
        int newPos = collapsed.pos();
        String newCigar = collapsed.cigar();

        if(isTxMatchPrimary(primaryRes))
        {
            OverhangGate.Result extended = mOverhangGate.extendSoftClips(
                    primaryRes.finalChromosome(), newPos, newCigar, readBases);
            newPos = extended.pos();
            newCigar = extended.cigar();
        }

        if(newPos != primaryRes.finalPos() || !newCigar.equals(primaryRes.finalCigar()))
        {
            TARS_LOGGER.trace(
                    "overhang-gate primary {}: {} -> {}", primary.getReadName(), primaryRes.finalCigar(), newCigar);
            primaryRes = primaryRes.withRevisedPrimary(newPos, newCigar, primaryRes.updatedMapQuality(), "overhang-gated");
            resolved[primaryIdx] = primaryRes;
        }

        // extend terminal softclips on tx-match XA alts too; genomic alts keep bwa's clip
        extendTxMatchAlts(primaryRes, primary);
    }

    private static byte[] reverseComplement(final byte[] bases)
    {
        byte[] reversed = Arrays.copyOf(bases, bases.length);
        Nucleotides.reverseComplementBasesInPlace(reversed, 0, reversed.length);
        return reversed;
    }

    // True when the chosen primary's placement was lifted from a transcript contig; gates the softclip extension.
    private static boolean isTxMatchPrimary(final LiftedRecord primaryRes)
    {
        LiftedAlignment primary = primaryRes.primaryAlignment();
        return primary != null && primary.FromTxContig;
    }

    // Extends a surviving terminal softclip on each tx-match XA alt; genomic alts and the primary are skipped.
    private void extendTxMatchAlts(final LiftedRecord primaryRes, final SAMRecord primary)
    {
        List<LiftedAlignment> alignments = primaryRes.liftedAlignments();
        int primaryIndex = primaryRes.primaryIndex();

        byte[] forwardBases = primary.getReadBases();
        if(forwardBases == null || forwardBases.length == 0)
        {
            return;
        }

        boolean recordForward = !primary.getReadNegativeStrandFlag();
        byte[] reverseBases = null;

        for(int i = 0; i < alignments.size(); ++i)
        {
            LiftedAlignment alt = alignments.get(i);
            if(i == primaryIndex || alt.Dropped || !alt.FromTxContig)
            {
                continue;
            }
            if(alt.LiftedCigar == null || alt.LiftedCigar.indexOf('S') < 0)
            {
                continue;
            }

            byte[] bases;
            if(alt.ForwardStrand == recordForward)
            {
                bases = forwardBases;
            }
            else
            {
                if(reverseBases == null)
                {
                    reverseBases = reverseComplement(forwardBases);
                }
                bases = reverseBases;
            }

            OverhangGate.Result extended = mOverhangGate.extendSoftClips(
                    alt.LiftedChromosome, alt.LiftedPos, alt.LiftedCigar, bases);
            if(extended.pos() != alt.LiftedPos || !extended.cigar().equals(alt.LiftedCigar))
            {
                alignments.set(i, alt.withLiftedCigar(extended.pos(), extended.cigar()));
            }
        }
    }
}
