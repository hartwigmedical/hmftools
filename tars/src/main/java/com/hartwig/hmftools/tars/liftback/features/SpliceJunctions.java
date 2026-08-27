package com.hartwig.hmftools.tars.liftback.features;

import java.nio.charset.StandardCharsets;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Rates a candidate intron against the annotated junctions from the sidecar and, failing that, the reference splice
// motif at its flanks.
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
        if(mRefGenome == null)
        {
            return Tier.NONE;
        }
        byte[] donor = refBases(
                candidateIntron.Chromosome, candidateIntron.start(), candidateIntron.start() + 1);
        byte[] acceptor = refBases(
                candidateIntron.Chromosome, candidateIntron.end() - 1, candidateIntron.end());
        return motifTier(donor, acceptor);
    }

    public int strand(final ChrBaseRegion intron)
    {
        int annotatedStrand = mEnsemblAnnotationIndex.junctionStrand(intron);
        if(annotatedStrand != 0 || mRefGenome == null)
        {
            return annotatedStrand;
        }

        byte[] donor = refBases(intron.Chromosome, intron.start(), intron.start() + 1);
        byte[] acceptor = refBases(intron.Chromosome, intron.end() - 1, intron.end());
        return motifStrand(donor, acceptor);
    }

    // 0 unless every N gap in the cigar agrees on a strand
    public int spliceStrand(final String chromosome, final int start, final String cigar)
    {
        int referencePosition = start;
        int strand = 0;
        for(CigarElement element : CigarUtils.cigarElementsFromStr(cigar))
        {
            if(element.getOperator() == CigarOperator.N)
            {
                int intronStrand = strand(new ChrBaseRegion(
                        chromosome, referencePosition, referencePosition + element.getLength() - 1));
                if(intronStrand == 0 || strand != 0 && strand != intronStrand)
                {
                    return 0;
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

    static int motifStrand(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2 || acceptorBases == null || acceptorBases.length != 2)
        {
            return 0;
        }

        String donor = upperCase(donorBases);
        String acceptor = upperCase(acceptorBases);
        if(donor.equals("GT") && acceptor.equals("AG")
                || donor.equals("GC") && acceptor.equals("AG")
                || donor.equals("AT") && acceptor.equals("AC"))
        {
            return 1;
        }
        if(donor.equals("CT") && acceptor.equals("AC")
                || donor.equals("CT") && acceptor.equals("GC")
                || donor.equals("GT") && acceptor.equals("AT"))
        {
            return -1;
        }
        return 0;
    }

    // Donor/acceptor 2-base flanks: GT-AG canonical (~99% of sites), GC-AG and AT-AC semi-canonical. The strand is
    // unknown at scan time, so each motif's reverse complement is accepted too.
    static Tier motifTier(final byte[] donorBases, final byte[] acceptorBases)
    {
        if(donorBases == null || donorBases.length != 2 || acceptorBases == null || acceptorBases.length != 2)
        {
            return Tier.NONE;
        }

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

    private byte[] refBases(final String chromosome, final int posStart, final int posEnd)
    {
        return BwaScoring.refWindow(mRefGenome, chromosome, posStart, posEnd);
    }
}
