package com.hartwig.hmftools.tars.liftback.merge;

import static com.hartwig.hmftools.tars.liftback.GenomicAlignmentScorer.referenceBases;

import java.nio.charset.StandardCharsets;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class SpliceJunctions
{
    // Splice-motif strength of a candidate junction, weakest to strongest. Declaration order is meaningful: candidates
    // are ranked with compareTo.
    public enum Tier
    {
        NONE,            // neither donor nor acceptor matches a known splice motif
        SEMI_CANONICAL,  // GC-AG / AT-AC (and reverse-complement equivalents)
        CANONICAL,       // GT-AG (~99% of splice sites)
        ANNOTATED        // matches an annotated junction from the sidecar
    }

    private final EnsemblAnnotationIndex mEnsemblAnnotationIndex;
    private final RefGenomeInterface mRefGenome;

    public SpliceJunctions(final EnsemblAnnotationIndex annotationIndex, final RefGenomeInterface refGenome)
    {
        mEnsemblAnnotationIndex = annotationIndex;
        mRefGenome = refGenome;
    }

    public Tier tier(final ChrBaseRegion candidateIntron)
    {
        if(mEnsemblAnnotationIndex.containsJunction(candidateIntron))
        {
            return Tier.ANNOTATED;
        }
        byte[] donor = referenceBases(mRefGenome, candidateIntron.Chromosome, candidateIntron.start(), candidateIntron.start() + 1);
        byte[] acceptor = referenceBases(mRefGenome, candidateIntron.Chromosome, candidateIntron.end() - 1, candidateIntron.end());
        return motifTier(donor, acceptor);
    }

    @Nullable
    public Orientation strand(final ChrBaseRegion intron)
    {
        Orientation annotatedStrand = mEnsemblAnnotationIndex.junctionStrand(intron);
        if(annotatedStrand != null)
        {
            return annotatedStrand;
        }

        byte[] donor = referenceBases(mRefGenome, intron.Chromosome, intron.start(), intron.start() + 1);
        byte[] acceptor = referenceBases(mRefGenome, intron.Chromosome, intron.end() - 1, intron.end());
        return motifStrand(donor, acceptor);
    }

    // 0 unless every N gap in the cigar agrees on a strand
    @Nullable
    public Orientation spliceStrand(final String chromosome, final int start, final String cigar)
    {
        int referencePosition = start;
        Orientation strand = null;
        for(CigarElement element : CigarUtils.cigarElementsFromStr(cigar))
        {
            if(element.getOperator() == CigarOperator.N)
            {
                Orientation intronStrand = strand(new ChrBaseRegion(
                        chromosome, referencePosition, referencePosition + element.getLength() - 1));
                if(intronStrand == null || strand != null && strand != intronStrand)
                {
                    return null;
                }
                strand = intronStrand;
            }
            if(element.getOperator().consumesReferenceBases())
            {
                referencePosition += element.getLength();
            }
        }
        return strand;
    }

    @Nullable
    static Orientation motifStrand(final byte[] donorBases, final byte[] acceptorBases)
    {
        String donor = upperCase(donorBases);
        String acceptor = upperCase(acceptorBases);
        if(donor.equals("GT") && acceptor.equals("AG")
                || donor.equals("GC") && acceptor.equals("AG")
                || donor.equals("AT") && acceptor.equals("AC"))
        {
            return Orientation.FORWARD;
        }
        if(donor.equals("CT") && acceptor.equals("AC")
                || donor.equals("CT") && acceptor.equals("GC")
                || donor.equals("GT") && acceptor.equals("AT"))
        {
            return Orientation.REVERSE;
        }
        return null;
    }

    // Donor/acceptor 2-base flanks: GT-AG canonical (~99% of sites), GC-AG and AT-AC semi-canonical. The strand is
    // unknown at scan time, so each motif's reverse complement is accepted too.
    static Tier motifTier(final byte[] donorBases, final byte[] acceptorBases)
    {
        String donor = upperCase(donorBases);
        String acceptor = upperCase(acceptorBases);

        if(donor.equals("GT") && acceptor.equals("AG") || donor.equals("CT") && acceptor.equals("AC"))
        {
            return Tier.CANONICAL;
        }
        if(donor.equals("GC") && acceptor.equals("AG") || donor.equals("CT") && acceptor.equals("GC"))
        {
            return Tier.SEMI_CANONICAL;
        }
        if(donor.equals("AT") && acceptor.equals("AC") || donor.equals("GT") && acceptor.equals("AT"))
        {
            return Tier.SEMI_CANONICAL;
        }
        return Tier.NONE;
    }

    // the reference is soft-masked, so repeat-region bases arrive lower case
    private static String upperCase(final byte[] bases)
    {
        return new String(bases, StandardCharsets.US_ASCII).toUpperCase();
    }
}
