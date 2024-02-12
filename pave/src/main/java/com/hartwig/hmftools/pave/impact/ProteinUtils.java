package com.hartwig.hmftools.pave.impact;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.aminoAcidFromBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pave.VariantData;

public final class ProteinUtils
{
    public static ProteinContext determineContext(
            final VariantData variant, final CodingContext cc, final TranscriptData transData, final RefGenomeInterface refGenome)
    {
        ProteinContext pc = new ProteinContext();

        boolean posStrand = transData.posStrand();

        ExonData exon = transData.exons().stream().filter(x -> x.Rank == cc.ExonRank).findFirst().orElse(null);

        if(exon == null)
        {
            PV_LOGGER.error("var({}) trans({}) invalid coding exonRank({})", variant, transData.TransName, cc.ExonRank);
            return pc;
        }

        // MNVs will have the coding range set exactly over their ref/alt bases (minus non-coding bases), whereas
        // DELs will cover the first ref base either side of the DEL, and inserts will be just to upstream ref base

        // both coding base and amino acid position start at 1, so eg coding base of 1-3 = amino acid 1, 4-6 = 2 etc
        // this CodonIndex may precede the first alt-base for INDELs for the reason described above
        pc.CodonIndex = (cc.CodingBase - 1) / 3 + 1;

        String refCodingBases = refGenome.getBaseString(variant.Chromosome, cc.CodingPositionRange[SE_START], cc.CodingPositionRange[SE_END]);
        pc.RefCodonBases = refCodingBases;

        pc.RefCodonsRanges.add(new int[] { cc.CodingPositionRange[SE_START], cc.CodingPositionRange[SE_END] });

        int upstreamOpenCodonBases = cc.CodingBase > 1 ? getOpenCodonBases(cc.UpstreamPhase) : 0;

        if(upstreamOpenCodonBases > 0)
        {
            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = !posStrand;

            String upstreamBases = getExtraBases(
                    transData, refGenome, variant.Chromosome, exon, upstreamOpenCodonBases, searchUp, pc.RefCodonsRanges);

            if(posStrand)
            {
                pc.RefCodonBases = upstreamBases + pc.RefCodonBases;
            }
            else
            {
                pc.RefCodonBases += upstreamBases;
            }
        }

        int downstreamMod = pc.RefCodonBases.length() % 3;
        int downstreamOpenCodonBases = downstreamMod == 0 ? 0 : 3 - downstreamMod;

        if(variant.isInsert())
            downstreamOpenCodonBases += 3; // get the ref codon past the insert as well

        if(downstreamOpenCodonBases > 0)
        {
            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = posStrand;

            String downstreamBases = getExtraBases(
                    transData, refGenome, variant.Chromosome, exon, downstreamOpenCodonBases, searchUp, pc.RefCodonsRanges);

            if(posStrand)
            {
                pc.RefCodonBases += downstreamBases;
            }
            else
            {
                pc.RefCodonBases = downstreamBases + pc.RefCodonBases;
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

            // MNV example: var GCT > AAT, with ref: ACG-CTA and alt: ACA-ATA
            // var length = 3, open codon bases = 1 (ie A), ref codon length = 6
            // working: alt = post-var (6 - 3 - 1) AC + alt AAT + upstream pre-var (6 - 1) A

            // DEL example: var TACA > T, with ref: CGT-ACA-ACA and alt: CGT-...-ACA
            // var length = 4, open codon bases = 2 (ie CA), ref codon length = 9, override upstream bases = 3
            // working: alt = post-var (9 - 4 - 3) CG + alt T (1st ref base) + upstream pre-var (9 - 3) ACA
            // but if the delete on the neg strand finishes at the exon boundary or coding region, it will already have been
            // truncated 1 base

            int upstreamRefBases = variant.isDeletion() && cc.CodingPositionRange[SE_END] == variant.EndPosition
                    ? upstreamOpenCodonBases + 1 : upstreamOpenCodonBases;

            if(variant.isInsert())
            {
                // adjustments for inserts since they have the ref stuck on the upper side
                alt = alt.substring(1) + refCodingBases;
            }
            else if(variant.isDeletion())
            {
                alt = !alt.isEmpty() ? alt.substring(0, 1) : "";
            }

            if(upstreamRefBases > 0 && codonBaseLength >= varLength + upstreamRefBases)
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength - upstreamRefBases) + alt
                        + pc.RefCodonBases.substring(codonBaseLength - upstreamRefBases);
            }
            else
            {
                pc.AltCodonBases = pc.RefCodonBases.substring(0, codonBaseLength - varLength) + alt;
            }
        }

        if(posStrand)
        {
            pc.RefAminoAcids = aminoAcidFromBases(pc.RefCodonBases);
        }
        else
        {
            pc.RefAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.RefCodonBases));
        }

        pc.AltCodonBasesComplete = pc.AltCodonBases;

        // fill in any incomplete alt AAs due to a frameshift
        int downstreamAltMod = pc.AltCodonBases.length() % 3;
        int downstreamAltOpenCodonBases = downstreamAltMod == 0 ? 0 : 3 - downstreamAltMod;

        if(downstreamAltOpenCodonBases > 0)
        {
            int downstreamStartPos = posStrand ? pc.refCodingBaseEnd() : pc.refCodingBaseStart();

            // get bases upstream to complete the upstream part of the codon
            boolean searchUp = posStrand;

            String downstreamBases = getExtraBases(
                    transData, refGenome, variant.Chromosome, exon, downstreamStartPos, downstreamAltOpenCodonBases, searchUp);

            if(posStrand)
            {
                pc.AltCodonBasesComplete += downstreamBases;
            }
            else
            {
                pc.AltCodonBasesComplete = downstreamBases + pc.AltCodonBasesComplete;
            }
        }

        if(!pc.validAltCodon())
        {
            pc.NetRefAminoAcids = pc.RefAminoAcids;
            pc.NetCodonIndexRange[SE_START] = pc.CodonIndex;
            pc.NetCodonIndexRange[SE_END] = pc.CodonIndex + pc.RefAminoAcids.length() - 1;
            return pc;
        }

        if(posStrand)
        {
            pc.AltAminoAcids = aminoAcidFromBases(pc.AltCodonBasesComplete);
        }
        else
        {
            pc.AltAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.AltCodonBasesComplete));
        }

        // extend the range of up and downstream codons for INDELs to help determine the net differences in amino acids
        extendIndelCodonBases(variant, cc, pc, transData, exon, refGenome, downstreamAltOpenCodonBases);

        trimAminoAcids(variant, pc);

        // check for a duplication of AAs
        checkInsertDuplication(variant, cc, pc, transData, exon, refGenome);

        return pc;
    }

    public static void extendIndelCodonBases(
            final VariantData variant, final CodingContext cc, final ProteinContext pc,
            final TranscriptData transData, final ExonData exon, final RefGenomeInterface refGenome, int downstreamAltOpenCodonBases)
    {
        if(!variant.isIndel())
            return;

        boolean posStrand = transData.posStrand();

        // for indels contained with a codon (ie phase 1 or 2), extract the first additional upstream bases
        if(variant.isInsert() && (cc.UpstreamPhase == PHASE_1 || cc.UpstreamPhase == PHASE_2) && pc.CodonIndex > 1)
        {
            boolean firstMatchesLast = pc.RefAminoAcids.charAt(0) == pc.AltAminoAcids.charAt(pc.AltAminoAcids.length() - 1);

            boolean firstSeqMatchesLast = pc.RefAminoAcids.length() >= 2 && pc.AltAminoAcids.length() >= 2
                    && pc.RefAminoAcids.substring(0, 2).equals(pc.AltAminoAcids.substring(pc.AltAminoAcids.length() - 2));

            if(firstMatchesLast || firstSeqMatchesLast)
            {
                boolean searchUp = !posStrand;

                String upstreamBases = getExtraBases(
                        transData, refGenome, variant.Chromosome, exon, 3, searchUp, pc.RefCodonsRanges);

                if(upstreamBases != null)
                {
                    --pc.CodonIndex;
                    pc.ExtraUpstreamCodon = true;

                    if(posStrand)
                    {
                        pc.RefCodonBases = upstreamBases + pc.RefCodonBases;
                        pc.AltCodonBases = upstreamBases + pc.AltCodonBases;
                        pc.AltCodonBasesComplete = upstreamBases + pc.AltCodonBasesComplete;

                        pc.RefAminoAcids = aminoAcidFromBases(pc.RefCodonBases);
                        pc.AltAminoAcids = aminoAcidFromBases(pc.AltCodonBasesComplete);
                    }
                    else
                    {
                        pc.RefCodonBases += upstreamBases;
                        pc.AltCodonBases += upstreamBases;
                        pc.AltCodonBasesComplete += upstreamBases;

                        pc.RefAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.RefCodonBases));
                        pc.AltAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.AltCodonBasesComplete));
                    }
                }
            }
        }

        if(!pc.RefAminoAcids.equals(pc.AltAminoAcids) || downstreamAltOpenCodonBases == 0)
            return;

        // from last ref codon base get 3 more each time
        boolean searchUp = posStrand;
        int currentPos = posStrand ? pc.refCodingBaseEnd() : pc.refCodingBaseStart();
        int extraBases = 3 + downstreamAltOpenCodonBases;

        while(true)
        {
            String downstreamBases = getExtraBases(
                    transData, refGenome, variant.Chromosome, exon, currentPos, extraBases, searchUp);

            if(downstreamBases == null)
                break;

            boolean codingRegionEnded = downstreamBases.length() < extraBases || downstreamBases.length() < downstreamAltOpenCodonBases;

            if(codingRegionEnded)
            {
                // occurs if the coding region end is reached
                PV_LOGGER.debug("var({}) invalid extraBases({}) downstreamBases({}) downstreamAltOpenCodonBases({})",
                        variant, extraBases, downstreamBases, downstreamAltOpenCodonBases);
            }

            if(posStrand)
            {
                // trim of the extra base(s) that the frame-shifted alt needs
                String refDownstreamBases = downstreamBases.substring(0, downstreamBases.length() - downstreamAltOpenCodonBases);
                pc.RefCodonBasesExtended = pc.RefCodonBases + refDownstreamBases;
                pc.AltCodonBasesComplete = pc.AltCodonBases + downstreamBases;
                pc.RefAminoAcids = aminoAcidFromBases(pc.RefCodonBasesExtended);
                pc.AltAminoAcids = aminoAcidFromBases(pc.AltCodonBasesComplete);
            }
            else
            {
                String refDownstreamBases = downstreamBases.substring(downstreamAltOpenCodonBases);
                pc.RefCodonBasesExtended = refDownstreamBases + pc.RefCodonBases;
                pc.AltCodonBasesComplete = downstreamBases + pc.AltCodonBases;
                pc.RefAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.RefCodonBasesExtended));
                pc.AltAminoAcids = aminoAcidFromBases(reverseComplementBases(pc.AltCodonBasesComplete));
            }

            if(!pc.RefAminoAcids.equals(pc.AltAminoAcids))
                break;

            if(codingRegionEnded)
                break;

            extraBases += 3;
        }
    }

    private static void trimAminoAcids(final VariantData variant, final ProteinContext proteinContext)
    {
        // trim matching amino acids from the start and end to align any differences to the most 3'UTR end in protein space
        if(proteinContext.RefAminoAcids.equals(proteinContext.AltAminoAcids))
            return; // ignore synonymous

        // determine features of the net impact on amino acids:
        // MNVs: strip off any synonymous AAs from either end
        // DEL: strip off the initial AA if it doesn't overlap a deleted base
        // INS:
        boolean canTrimStart = variant.isBaseChange() || variant.isDeletion() || variant.isInsert();

        boolean repeatStartRemoval = canTrimStart && variant.isIndel();

        boolean repeatEndRemoval = variant.isInsert();

        trimAminoAcids(proteinContext, canTrimStart, repeatStartRemoval, repeatEndRemoval);
    }

    public static void trimAminoAcids(
            final ProteinContext proteinContext, boolean canTrimStart, boolean repeatStartRemoval, boolean repeatEndRemoval)
    {
        int aaIndex = proteinContext.CodonIndex;
        String netRefAminoAcids = proteinContext.RefAminoAcids;
        String netAltAminoAcids = proteinContext.AltAminoAcids;

        // strip off the initial AA if it doesn't relate to any of the altered bases, as is the case for a DEL with the reference
        // is the last base of a codon
        if(canTrimStart)
        {
            // strip out any matching AA from the start if the ref has more than 1
            while(true)
            {
                if(netRefAminoAcids.length() > 1 && !netRefAminoAcids.isEmpty() && !netAltAminoAcids.isEmpty()
                && netRefAminoAcids.charAt(0) == netAltAminoAcids.charAt(0))
                {
                    netRefAminoAcids = netRefAminoAcids.substring(1);
                    netAltAminoAcids = netAltAminoAcids.substring(1);
                    aaIndex++;

                    if(!repeatStartRemoval) // only the first
                        break;
                }
                else
                {
                    break;
                }
            }
        }

        // and strip off a matching AA from the end if there is one
        while(true)
        {
            if(!netRefAminoAcids.isEmpty() && !netAltAminoAcids.isEmpty()
            && netRefAminoAcids.charAt(netRefAminoAcids.length() - 1) == netAltAminoAcids.charAt(netAltAminoAcids.length() - 1))
            {
                netRefAminoAcids = netRefAminoAcids.substring(0, netRefAminoAcids.length() - 1);
                netAltAminoAcids = netAltAminoAcids.substring(0, netAltAminoAcids.length() - 1);

                if(!repeatEndRemoval)
                    break;
            }
            else
            {
                break;
            }
        }

        proteinContext.NetRefAminoAcids = netRefAminoAcids;
        proteinContext.NetAltAminoAcids = netAltAminoAcids;
        proteinContext.NetCodonIndexRange[SE_START] = aaIndex;
        proteinContext.NetCodonIndexRange[SE_END] = aaIndex + max(netRefAminoAcids.length() - 1, 0);
    }

    public static void checkInsertDuplication(
            final VariantData variant, final CodingContext cc, final ProteinContext pc,
            final TranscriptData transData, final ExonData exon, final RefGenomeInterface refGenome)
    {
        if(!variant.isInsert() || cc.IsFrameShift)
            return;

        if(pc.NetAltAminoAcids.isEmpty() || pc.RefAminoAcids.length() < 2)
            return;

        // first inserted amino acids matches end ref AA
        if(pc.NetAltAminoAcids.charAt(0) != pc.RefAminoAcids.charAt(0) && pc.NetAltAminoAcids.charAt(0) != pc.RefAminoAcids.charAt(1))
            return;

        // a duplication of AAs means 1 extra copy, not multiple
        // duplication (single AA) p.Gln8dup
        // duplication (range) p.Gly4_Gln6dup - variety of AAs
        // dup range, single AA repeated: c.1014_1019dupTGCTGC	p.Ala339_Ala340dup

        boolean posStrand = transData.posStrand();

        // check for a duplication of AAs

        // check for a single AA repeated
        final String netAminoAcids = pc.NetAltAminoAcids;

        if(pc.NetAltAminoAcids.length() == 1 && pc.NetRefAminoAcids.isEmpty())
        {
            if(pc.NetAltAminoAcids.charAt(0) == pc.RefAminoAcids.charAt(0))
            {
                pc.IsDuplication = true;
                pc.NetCodonIndexRange[SE_END] = pc.NetCodonIndexRange[SE_START] = pc.CodonIndex;
                return;
            }
            else if(pc.NetAltAminoAcids.charAt(0) == pc.RefAminoAcids.charAt(1))
            {
                pc.IsDuplication = true;
                pc.NetCodonIndexRange[SE_END] = pc.NetCodonIndexRange[SE_START];
                return;
            }
        }

        // otherwise extend downstream the length of the inserted AAs and test for a match
        int extraBases = (netAminoAcids.length() - 1) * 3;
        boolean searchUp = posStrand;
        int currentPos = posStrand ? pc.refCodingBaseEnd() : pc.refCodingBaseStart();

        String downstreamBases = getExtraBases(
                transData, refGenome, variant.Chromosome, exon, currentPos, extraBases, searchUp);

        if(downstreamBases == null)
            return;

        String extraAminoAcids = posStrand
                ? aminoAcidFromBases(downstreamBases)
                : aminoAcidFromBases(reverseComplementBases(downstreamBases));

        String extendedRefAminoAcids = pc.RefAminoAcids.charAt(1) + extraAminoAcids;

        if(extendedRefAminoAcids.equals(netAminoAcids))
        {
            pc.IsDuplication = true;
            pc.NetCodonIndexRange[SE_END] = pc.NetCodonIndexRange[SE_START] + extendedRefAminoAcids.length() - 1;
        }
    }

    public static int getOpenCodonBases(int phase)
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
        // convenience method which doesn't require a list of coding base ranges, nor updates them
        final List<int[]> codingBaseRanges = Lists.newArrayList(new int[] {startPos, startPos});
        return getExtraBases(transData, refGenome, chromosome, currentExon, requiredBases, searchUp, codingBaseRanges);
    }

    public static String getExtraBases(
            final TranscriptData transData, final RefGenomeInterface refGenome, final String chromosome,
            final ExonData refExon, int requiredBases, boolean searchUp, final List<int[]> codingBaseRanges)
    {
        if(codingBaseRanges.isEmpty())
            return null;

        String extraBases = "";

        int[] currentRange = searchUp ? codingBaseRanges.get(codingBaseRanges.size() - 1) : codingBaseRanges.get(0);

        int startPos = searchUp ? currentRange[SE_END] : currentRange[SE_START];

        final ExonData currentExon = positionWithin(startPos, refExon.Start, refExon.End) ?
                refExon : transData.exons().stream().filter(x -> positionWithin(startPos, x.Start, x.End)).findFirst().orElse(null);

        if(currentExon == null)
            return null;

        if(searchUp)
        {
            int currentExonBases = min(requiredBases, min(currentExon.End, transData.CodingEnd) - startPos);

            if(currentExonBases > 0)
            {
                int upperPos = startPos + 1 + (currentExonBases - 1);
                currentRange[SE_END] = upperPos;

                extraBases = refGenome.getBaseString(chromosome, startPos + 1, upperPos);
                requiredBases -= currentExonBases;
            }

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank + 1 : currentExon.Rank - 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                if(nextExon == null || nextExon.Start > transData.CodingEnd)
                    return extraBases;

                int lowerPos = nextExon.Start;
                int upperPos = nextExon.Start + (requiredBases - 1);
                int[] newRange = { lowerPos, upperPos };
                codingBaseRanges.add(newRange);

                String nextExonBases = refGenome.getBaseString(chromosome, lowerPos, upperPos);
                extraBases += nextExonBases;
            }
        }
        else
        {
            // search in lower positions
            int currentExonBases = min(requiredBases, startPos - max(currentExon.Start, transData.CodingStart));

            if(currentExonBases > 0)
            {
                int lowerPos = startPos - currentExonBases;
                currentRange[SE_START] = lowerPos;

                extraBases = refGenome.getBaseString(chromosome, lowerPos, startPos - 1);
                requiredBases -= currentExonBases;
            }

            if(requiredBases > 0)
            {
                int nextExonRank = transData.posStrand() ? currentExon.Rank - 1 : currentExon.Rank + 1;

                ExonData nextExon = transData.exons().stream().filter(x -> x.Rank == nextExonRank).findFirst().orElse(null);

                if(nextExon == null || nextExon.End < transData.CodingStart)
                    return extraBases;

                int lowerPos = nextExon.End - (requiredBases - 1);
                int upperPos = nextExon.End;
                int[] newRange = { lowerPos, upperPos };
                codingBaseRanges.add(0, newRange);

                String nextExonBases = refGenome.getBaseString(chromosome, lowerPos, upperPos);
                extraBases = nextExonBases + extraBases;
            }
        }

        return extraBases;
    }
}
