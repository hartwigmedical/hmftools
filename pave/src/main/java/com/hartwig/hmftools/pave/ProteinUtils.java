package com.hartwig.hmftools.pave;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public final class ProteinUtils
{
    public static ProteinContext determineContext(
            final VariantData variant, final CodingContext cc, final TranscriptData transData, final RefGenomeInterface refGenome)
    {
        ProteinContext pc = new ProteinContext();

        // both coding base and amino acid position start at 1
        // coding base of 1-3 = amino acid 1, 4-6 = 2 etc
        pc.CodonIndex = (cc.CodingBase - 1) / 3 + 1;

        boolean posStrand = transData.posStrand();

        int upstreamStartPos = posStrand ? cc.CodingPositionRange[SE_START] : cc.CodingPositionRange[SE_END];

        ExonData exon = transData.exons().stream().filter(x -> x.Rank == cc.ExonRank).findFirst().orElse(null);

        if(exon == null)
        {
            PV_LOGGER.error("var({}) trans({}) invalid coding exonRank({})", variant, transData.TransName, cc.ExonRank);
            return pc;
        }

        String refCodingBases = refGenome.getBaseString(variant.Chromosome, cc.CodingPositionRange[SE_START], cc.CodingPositionRange[SE_END]);
        pc.RefCodonBases = refCodingBases;

        pc.RefCodonsRange[SE_START] = cc.CodingPositionRange[SE_START];
        pc.RefCodonsRange[SE_END] = cc.CodingPositionRange[SE_END];

        int upstreamOpenCodonBases = getOpenCodonBases(cc.UpstreamPhase);

        if(upstreamOpenCodonBases > 0)
        {
            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = !posStrand;
            String upstreamBases = getExtraBases(transData, refGenome, variant.Chromosome, exon, upstreamStartPos, upstreamOpenCodonBases, searchUp);

            if(posStrand)
            {
                pc.RefCodonBases = upstreamBases + pc.RefCodonBases;
                pc.RefCodonsRange[SE_START] -= upstreamOpenCodonBases;
            }
            else
            {
                pc.RefCodonBases += upstreamBases;
                pc.RefCodonsRange[SE_END] += upstreamOpenCodonBases;
            }
        }

        int downstreamMod = pc.RefCodonBases.length() % 3;
        int downstreamOpenCodonBases = downstreamMod == 0 ? 0 : 3 - downstreamMod;

        if(downstreamOpenCodonBases > 0)
        {
            //int downstreamStartPos = posStrand ? variant.EndPosition : variant.Position;
            int downstreamStartPos = posStrand ? cc.CodingPositionRange[SE_END] : cc.CodingPositionRange[SE_START];

            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = posStrand;
            String downstreamBases = getExtraBases(transData, refGenome, variant.Chromosome, exon, downstreamStartPos, downstreamOpenCodonBases, searchUp);

            if(posStrand)
            {
                pc.RefCodonBases += downstreamBases;
                pc.RefCodonsRange[SE_END] += downstreamOpenCodonBases;
            }
            else
            {
                pc.RefCodonBases = downstreamBases + pc.RefCodonBases;
                pc.RefCodonsRange[SE_START] -= downstreamOpenCodonBases;
            }
        }

        // if codon(s) are incomplete then either a bug or an issue with the transcript definition
        if(!pc.validRefCodon())
            return pc;

        // only factor in coding bases in the mutation
        String ref = cc.codingRef(variant);
        String alt = cc.codingAlt(variant);
        int varLength = ref.length();

        // find the start of the codon in which the first base of the variant is in, and likewise the end of the last codon

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
        }
        else
        {
            int codonBaseLength = pc.RefCodonBases.length();

            // adjustments for inserts since they have the ref stuck on the upper side
            if(variant.isInsert())
            {
                alt = alt.substring(1) + refCodingBases;
            }

            if(upstreamOpenCodonBases > 0)
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength - upstreamOpenCodonBases) + alt
                        + pc.RefCodonBases.substring(codonBaseLength - upstreamOpenCodonBases);
            }
            else
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength) + alt;
            }
        }

        // an INDEL causing a frameshift can exit at this point since the novel AAs do not need to be recorded
        // (and may gone on until the end of the transcript)
        if(variant.isIndel())
        {
            int adjustedCodingBases = variant.isInsert() ? variant.baseDiff() : ref.length() - alt.length();

            if((adjustedCodingBases % 3) != 0)
            {
                cc.IsFrameShift = true;
                return pc;
            }
        }

        if(posStrand)
        {
            pc.RefAminoAcids = Codons.aminoAcidFromBases(pc.RefCodonBases);
            pc.AltAminoAcids = Codons.aminoAcidFromBases(pc.AltCodonBases);
        }
        else
        {
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

    public static String getExtraBases(
            final TranscriptData transData, final RefGenomeInterface refGenome, final String chromosome,
            final ExonData currentExon, int startPos, int requiredBases, boolean searchUp)
    {
        String extraBases = "";

        if(searchUp)
        {
            int currentExonBases = min(requiredBases, currentExon.End - startPos);

            if(currentExonBases > 0)
            {
                extraBases = refGenome.getBaseString(chromosome, startPos + 1, startPos + 1 + (currentExonBases - 1));
                requiredBases -= currentExonBases;
            }

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank + 1 : currentExon.Rank - 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                if(nextExon == null)
                    return extraBases;

                String nextExonBases = refGenome.getBaseString(chromosome, nextExon.Start, nextExon.Start + (requiredBases - 1));
                extraBases += nextExonBases;
            }
        }
        else
        {
            // search in lower positions
            int currentExonBases = min(requiredBases, startPos - currentExon.Start);

            if(currentExonBases > 0)
            {
                extraBases = refGenome.getBaseString(chromosome, startPos - currentExonBases, startPos - 1);
                requiredBases -= currentExonBases;
            }

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank - 1 : currentExon.Rank + 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                if(nextExon == null)
                    return extraBases;

                String nextExonBases = refGenome.getBaseString(chromosome, nextExon.End - (requiredBases - 1), nextExon.End);
                extraBases = nextExonBases + extraBases;
            }
        }

        return extraBases;
    }
}
