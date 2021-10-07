package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.CodingUtils.getExtraBases;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class ProteinContext
{
    public String RefCodonBases; // coding bases rounded expanded to cover whole codons
    public String AltCodonBases; // as above but with ref swapped for alt

    public int StartPosition; // first amino acid affected
    public String RefAminoAcids;
    public String AltAminoAcids;

    public ProteinContext()
    {
        RefCodonBases = "";
        AltCodonBases = "";

        StartPosition = 0;
        RefAminoAcids = "";
        AltAminoAcids = "";
    }

    public boolean hasCodingBases() { return !RefAminoAcids.isEmpty(); }

    public boolean hasProteinChange() { return !RefAminoAcids.equals(AltAminoAcids); }

    public String hgvsStr() { return "tbc"; }

    public static ProteinContext determineContext(
            final VariantData variant, final CodingContext cc, final TranscriptData transData, final RefGenomeInterface refGenome)
    {
        ProteinContext pc = new ProteinContext();

        // both coding base and amino acid position start at 1
        // coding base of 1-3 = amino acid 1, 4-6 = 2 etc
        pc.StartPosition = (cc.CodingBase - 1) / 3 + 1;

        boolean posStrand = transData.posStrand();
        int upstreamStartPos = posStrand ? variant.Position : variant.EndPosition;

        ExonData exon = transData.exons().stream().filter(x -> x.Rank == cc.ExonRank).findFirst().orElse(null);

        if(exon == null)
        {
            PV_LOGGER.error("var({}) trans({}) invalid coding exonRank({})", variant, transData.TransName, cc.ExonRank);
            return pc;
        }

        pc.RefCodonBases = refGenome.getBaseString(variant.Chromosome, cc.CodingPositionRange[SE_START], cc.CodingPositionRange[SE_END]);

        int upstreamOpenCodonBases = getOpenCodonBases(cc.UpstreamPhase);

        if(upstreamOpenCodonBases > 0)
        {
            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = !posStrand;
            String upstreamBases = getExtraBases(transData, refGenome, variant.Chromosome, exon, upstreamStartPos, upstreamOpenCodonBases, searchUp);

            if(posStrand)
                pc.RefCodonBases = upstreamBases + pc.RefCodonBases;
            else
                pc.RefCodonBases += upstreamBases;
        }

        int downstreamOpenCodonBases = getDownstreamOpenCodonBases(cc.UpstreamPhase, transData.Strand, variant.baseDiff(), variant.Alt);
        int downstreamStartPos = posStrand ? variant.EndPosition : variant.Position;

        if(downstreamOpenCodonBases > 0)
        {
            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = posStrand;
            String downstreamBases = getExtraBases(transData, refGenome, variant.Chromosome, exon, downstreamStartPos, downstreamOpenCodonBases, searchUp);

            if(posStrand)
                pc.RefCodonBases += downstreamBases;
            else
                pc.RefCodonBases = downstreamBases + pc.RefCodonBases;
        }

        // an INDEL causing a frameshift can exit at this point since the novel AAs do not need to be recorded
        // (and may gone on until the end of the transcript)
        if(variant.isIndel())
        {
            if(variant.isInsert() && (variant.baseDiff() % 3) != 3)
                return pc;

            // usually would be just a check of whether deleted bases % 3, but need to consider some being beyond the coding region
            if(variant.isDeletion() && (cc.DeletedCodingBases % 3) != 3)
                return pc;
        }

        // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon
        String alt = variant.Alt;
        int varLength = variant.Ref.length();

        if(posStrand)
        {
            if(upstreamOpenCodonBases > 0)
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, upstreamOpenCodonBases) + alt
                        + pc.RefCodonBases.substring(upstreamOpenCodonBases + varLength);
            }
            else
            {
                pc.AltCodonBases = alt + pc.RefCodonBases.substring(varLength);
            }

            pc.RefAminoAcids = Codons.aminoAcidFromBases(pc.RefCodonBases);
            pc.AltAminoAcids = Codons.aminoAcidFromBases(pc.AltCodonBases);
        }
        else
        {
            // TODO
            // int codingBaseLen = (int)(Math.ceil((upstreamOpenCodonBases + varLength) / 3.0)) * 3;
            int codonBaseLength = pc.RefCodonBases.length();

            if(upstreamOpenCodonBases > 0)
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength - upstreamOpenCodonBases) + alt
                        + pc.RefCodonBases.substring(codonBaseLength - upstreamOpenCodonBases);
            }
            else
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength) + alt;
            }

            pc.RefAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(pc.RefCodonBases));
            pc.AltAminoAcids = Codons.aminoAcidFromBases(Nucleotides.reverseStrandBases(pc.AltCodonBases));
        }

        return pc;
    }

    private static int getOpenCodonBases(int phase)
    {
        // determine the number of bases upstream to make a codon
        if(phase == PHASE_1)
            return 0;
        else if(phase == PHASE_2)
            return 1;
        else
            return 2;
    }

    public static int getDownstreamOpenCodonBases(int startPhase, int strand, int baseChange, final String alt)
    {
        // determine the phase at the base after the mutation - last base of an MNV/SNV, and next base for an INS or DEL
        int mutationTicks;

        if(strand == POS_STRAND)
        {
            mutationTicks = baseChange < 0 ? 1 : alt.length();
        }
        else
        {
            // need to go to the phase (prior to?) the ALT
            if(baseChange == 0)
                mutationTicks = alt.length();
            else if(baseChange > 0)
                mutationTicks = alt.length() + 1;
            else
                mutationTicks = 2;
        }

        int postMutationPhase = tickPhaseForward(startPhase, mutationTicks);

        if(postMutationPhase == PHASE_0)
            return 1;
        else if(postMutationPhase == PHASE_1)
            return 0;
        else
            return 2;
    }

}
