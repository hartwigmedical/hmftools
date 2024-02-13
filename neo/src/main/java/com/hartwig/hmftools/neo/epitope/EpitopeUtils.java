package com.hartwig.hmftools.neo.epitope;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.isStopCodon;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.codon.AminoAcidRna.AA_SELENOCYSTEINE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.codon.AminoAcidRna.STOP_SYMBOL;
import static com.hartwig.hmftools.common.codon.AminoAcidRna.convertDnaCodonToAminoAcid;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.neo.PeptideData;

import org.apache.commons.compress.utils.Lists;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class EpitopeUtils
{
    public static void setTranscriptContext(final NeoEpitope neData, final TranscriptData transData, int position, int stream)
    {
        // determine phasing, coding and region context
        boolean isUpstream = stream == FS_UP;

        if(position < transData.TransStart)
        {
            neData.ExonRank[stream] = transData.posStrand() ? 0 : transData.exons().size() - 1;
            neData.RegionType[stream] = transData.posStrand() ? UPSTREAM : DOWNSTREAM;
            neData.Phases[stream] = PHASE_NONE;
            return;
        }
        else if(position > transData.TransEnd)
        {
            neData.ExonRank[stream] = transData.posStrand() ? transData.exons().size() - 1 : 0;
            neData.RegionType[stream] = transData.posStrand() ? DOWNSTREAM : UPSTREAM;
            neData.Phases[stream] = PHASE_NONE;
            return;
        }

        for(ExonData exon : transData.exons())
        {
            if(position > exon.End)
                continue;

            if(position < exon.Start)
            {
                // intronic
                neData.RegionType[stream] = INTRONIC;

                // upstream pos-strand, before next exon then take the exon before's rank
                if(transData.Strand == POS_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.Rank - 1;
                else if(transData.Strand == NEG_STRAND && isUpstream)
                    neData.ExonRank[stream] = exon.Rank;
                else if(transData.Strand == POS_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.Rank;
                else if(transData.Strand == NEG_STRAND && !isUpstream)
                    neData.ExonRank[stream] = exon.Rank - 1;

                if(transData.Strand == POS_STRAND)
                    neData.Phases[stream] = exon.PhaseStart;
                else
                    neData.Phases[stream] = exon.PhaseEnd;
            }
            else if(positionWithin(position, exon.Start, exon.End))
            {
                neData.RegionType[stream] = EXONIC;
                neData.ExonRank[stream] = exon.Rank;
            }

            break;
        }
    }

    public static void setTranscriptCodingData(
            final NeoEpitope neData, final TranscriptData transData, int position, int insSeqLength, int stream)
    {
        if(transData.CodingStart != null)
        {
            if(positionWithin(position, transData.CodingStart, transData.CodingEnd))
                neData.CodingType[stream] = CODING;
            else if(transData.Strand == POS_STRAND && position < transData.CodingStart)
                neData.CodingType[stream] = UTR_5P;
            else if(transData.Strand == NEG_STRAND && position > transData.CodingEnd)
                neData.CodingType[stream] = UTR_5P;
            else
                neData.CodingType[stream] = UTR_3P;
        }
        else
        {
            neData.CodingType[stream] = NON_CODING;
            neData.Phases[stream] = -1;
        }

        if(neData.CodingType[stream] == CODING && neData.RegionType[stream] == EXONIC)
        {
            final CodingBaseData cbData = calcCodingBases(transData, position);
            neData.Phases[stream] = tickPhaseForward(cbData.Phase, insSeqLength);
        }
    }

    public static int getUpstreamOpenCodonBases(int startPhase, int strand, boolean isIndel)
    {
        int phase = startPhase;

        // return the number of open codon bases up to but not including the variant's position
        if(strand == POS_STRAND || !isIndel)
        {
            if(phase == PHASE_1)
                return 0;
            else if(phase == PHASE_2)
                return 1;
            else
                return isIndel ? 0 : 2; // INDEL leaves the mutation base unchanged, so if it ends a codon, it's not novel
        }
        else
        {
            if(phase == PHASE_1)
                return 1;
            else if(phase == PHASE_2)
                return 2;
            else
                return isIndel ? 0 : 2;
        }
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

    public static String getUpstreamCodingBases(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation, int requiredBases)
    {
        final CodingBaseExcerpt cbData = getUpstreamCodingBaseExcerpt(
                refGenome, transData, chromosome, nePosition, neOrientation, requiredBases);

        return cbData.Bases;
    }

    public static CodingBaseExcerpt getUpstreamCodingBaseExcerpt(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation, int requiredBases)
    {
        if(requiredBases <= 0 || transData.CodingStart == null)
            return null;

        if(!positionWithin(nePosition, transData.CodingStart, transData.CodingEnd))
            return null;

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        Cigar cigar = new Cigar();
        int lastMatchBase = -1;
        int startPos;
        int endPos;

        if(neOrientation == NEG_ORIENT)
        {
            startPos = 0;
            endPos = 0;

            // walk forwards through the exons, collecting up the required positions and bases

            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(nePosition > exon.End)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break; // no more coding bases

                int exonBaseStart = max(exon.Start, nePosition);
                int exonBaseEnd = min(transData.CodingEnd, exon.End);
                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonBaseCount)
                {
                    // take them all
                    baseStart = exonBaseStart;
                    baseEnd = exonBaseEnd;
                    requiredBases -= exonBaseCount;
                }
                else
                {
                    baseStart = exonBaseStart;
                    baseEnd = baseStart + requiredBases - 1;
                    requiredBases = 0;
                }

                baseString += refGenome.getBaseString(chromosome, baseStart, baseEnd);

                if(startPos == 0)
                    startPos = exonBaseStart;

                endPos = baseEnd;

                if(lastMatchBase > 0)
                    cigar.add(new CigarElement(baseStart - lastMatchBase - 1, CigarOperator.N));

                lastMatchBase = baseEnd;

                cigar.add(new CigarElement(baseEnd - baseStart + 1, CigarOperator.M));

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            startPos = 0;
            endPos = 0;

            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(nePosition < exon.Start)
                    continue;

                if(exon.End < transData.CodingStart)
                    break;

                int exonBaseEnd = min(exon.End, nePosition);
                int exonBaseStart = max(transData.CodingStart, exon.Start);
                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonBaseCount)
                {
                    // take them all
                    baseStart = exonBaseStart;
                    baseEnd = exonBaseEnd;
                    requiredBases -= exonBaseCount;
                }
                else
                {
                    baseEnd = exonBaseEnd;
                    baseStart = baseEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                baseString = refGenome.getBaseString(chromosome, baseStart, baseEnd) + baseString;

                if(endPos == 0)
                    endPos = exonBaseEnd;

                startPos = baseStart;

                if(lastMatchBase > 0)
                    cigar.add(new CigarElement(lastMatchBase - baseEnd -  - 1, CigarOperator.N));

                lastMatchBase = baseStart;

                cigar.add(new CigarElement(baseEnd - baseStart + 1, CigarOperator.M));

                if (requiredBases <= 0)
                    break;
            }

            final List<CigarElement> elements = cigar.getCigarElements();

            cigar = new Cigar();
            for(int i = elements.size() - 1; i >= 0; --i)
            {
                cigar.add(elements.get(i));
            }
        }

        return new CodingBaseExcerpt(baseString, startPos, endPos, cigar);
    }

    public static String getDownstreamCodingBases(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation, int requiredBases,
            boolean canStartInExon, boolean reqSpliceAcceptor, boolean reqAllBases)
    {
        final CodingBaseExcerpt cbData = getDownstreamCodingBaseExcerpt(
                refGenome, transData, chromosome, nePosition, neOrientation, requiredBases, canStartInExon, reqSpliceAcceptor, reqAllBases);

        return cbData.Bases;
    }

    public static CodingBaseExcerpt getDownstreamCodingBaseExcerpt(
            final RefGenomeInterface refGenome, final TranscriptData transData,
            final String chromosome, int nePosition, byte neOrientation, int requiredBases,
            boolean canStartInExon, boolean reqSpliceAcceptor, boolean reqAllBases)
    {
        if(requiredBases == 0)
            return null;

        int codingStart = transData.CodingStart != null && nePosition >= transData.CodingStart ? transData.CodingStart : transData.TransStart;
        int codingEnd = transData.CodingEnd != null && nePosition <= transData.CodingEnd ? transData.CodingEnd : transData.TransEnd;

        // if the position falls in the stop codon, take at least the last base
        if((transData.Strand == POS_STRAND && nePosition >= transData.TransEnd) || (transData.Strand == NEG_STRAND && nePosition <= transData.TransStart))
            return null;

        final List<ExonData> exonDataList = transData.exons();

        String baseString = "";

        Cigar cigar = new Cigar();
        int lastMatchBase = -1;
        int startPos;
        int endPos;

        if(neOrientation == NEG_ORIENT)
        {
            startPos = 0;
            endPos = 0;

            for(int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                // starts after the first exon, ie at the first splice acceptor
                if(reqSpliceAcceptor && exon.Rank == 1)
                    continue;

                if(nePosition > exon.End)
                    continue;

                if(positionWithin(nePosition, exon.Start, exon.End) && !canStartInExon)
                    continue; // will start at the next exon

                if(exon.Start > codingEnd && !reqAllBases)
                    break; // no more coding bases

                int exonBaseStart = max(nePosition, exon.Start);
                int exonBaseEnd = !reqAllBases ? min(codingEnd, exon.End) : exon.End;
                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonBaseCount || reqAllBases)
                {
                    // take them all
                    baseStart = exonBaseStart;
                    baseEnd = exonBaseEnd;

                    if(!reqAllBases)
                        requiredBases -= exonBaseCount;
                }
                else
                {
                    baseStart = exonBaseStart;
                    baseEnd = baseStart + requiredBases - 1;
                    requiredBases = 0;
                }

                baseString += refGenome.getBaseString(chromosome, baseStart, baseEnd);

                if(startPos == 0)
                    startPos = exonBaseStart;

                endPos = baseEnd;

                if(lastMatchBase > 0)
                    cigar.add(new CigarElement(baseStart - lastMatchBase - 1, CigarOperator.N));

                lastMatchBase = baseEnd;

                cigar.add(new CigarElement(baseEnd - baseStart + 1, CigarOperator.M));

                if (requiredBases <= 0 && !reqAllBases)
                    break;
            }
        }
        else
        {
            startPos = 0;
            endPos = 0;

            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(reqSpliceAcceptor && exon.Rank == 1)
                    continue;
                
                if(nePosition < exon.Start)
                    continue;

                if(codingEnd < exon.Start && !reqAllBases)
                    continue;

                if(positionWithin(nePosition, exon.Start, exon.End) && !canStartInExon)
                    continue;

                if(exon.End < codingStart && !reqAllBases)
                    break;

                int exonBaseEnd = min(nePosition, exon.End);

                int exonBaseStart = !reqAllBases ? max(codingStart, exon.Start) : exon.Start;
                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                int baseStart, baseEnd;

                if(requiredBases >= exonBaseCount || reqAllBases)
                {
                    // take them all
                    baseStart = exonBaseStart;
                    baseEnd = exonBaseEnd;

                    if(!reqAllBases)
                        requiredBases -= exonBaseCount;
                }
                else
                {
                    baseEnd = exonBaseEnd;
                    baseStart = baseEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                baseString = refGenome.getBaseString(chromosome, baseStart, baseEnd) + baseString;

                if(endPos == 0)
                    endPos = exonBaseEnd;

                startPos = baseStart;

                if(lastMatchBase > 0)
                    cigar.add(new CigarElement(lastMatchBase - baseEnd - 1, CigarOperator.N));

                lastMatchBase = baseStart;

                cigar.add(new CigarElement(baseEnd - baseStart + 1, CigarOperator.M));

                if (requiredBases <= 0 && !reqAllBases)
                    break;
            }

            final List<CigarElement> elements = cigar.getCigarElements();

            cigar = new Cigar();
            for(int i = elements.size() - 1; i >= 0; --i)
            {
                cigar.add(elements.get(i));
            }
        }

        return new CodingBaseExcerpt(baseString, startPos, endPos, cigar);
    }

    public static void adjustCodingBasesForStrand(final NeoEpitope neData)
    {
        // upstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // downstream strand 1, bases will be retrieved from left to right (lower to higher), no need for any conversion
        // upstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert
        // downstream strand -1, bases will be retrieved from left to right (lower to higher), need to reverse and convert

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(neData.strand(fs) == NEG_STRAND)
                neData.CodingBases[fs] = reverseComplementBases(neData.RawCodingBases[fs]);
            else
                neData.CodingBases[fs] = neData.RawCodingBases[fs];
        }
    }

    public static int calcNonMediatedDecayBases(final NeoEpitope neData)
    {
        // distance from (novel) stop codon to last splice acceptor
        if(!neData.NovelAcid.endsWith(STOP_SYMBOL))
            return -1;

        if(neData.NovelAcid.equals(STOP_SYMBOL)) // ignore stop-gained
            return -1;

        final TranscriptData transData = neData.TransData[FS_DOWN];
        final List<ExonData> exonDataList = transData.exons();
        int refPosition = neData.position(FS_DOWN);

        int exonicBaseCount = 0;

        if(neData.orientation(FS_DOWN) == NEG_ORIENT)
        {
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == exonDataList.size() - 1)
                    break;

                if (refPosition > exon.End)
                    continue;

                exonicBaseCount += exon.End - max(refPosition, exon.Start) + 1;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(i == 0)
                    break;

                if(refPosition < exon.Start)
                    continue;

                exonicBaseCount += min(refPosition, exon.End) - exon.Start + 1;
            }
        }

        int newCodingBases = neData.NovelAcid.length() * 3;

        if(exonicBaseCount >= newCodingBases)
            return exonicBaseCount - newCodingBases;

        return -1;
    }

    public static int calcStopCodonBases(final NeoEpitope neData)
    {
        // distance from the downstream position to the stop codon or the end of the transcript if there is none
        if(neData.NovelAcid.endsWith(STOP_SYMBOL))
            return neData.NovelCodonBases.length();

        final TranscriptData transData = neData.TransData[FS_DOWN];
        final List<ExonData> exonDataList = transData.exons();
        int refPosition = neData.position(FS_DOWN);

        int codingBaseCount = 0;

        if(neData.orientation(FS_DOWN) == NEG_ORIENT)
        {
            int codingStop = transData.CodingEnd != null ? transData.CodingEnd : transData.TransEnd;

            if(codingStop <= refPosition)
                return 0;

            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (refPosition > exon.End)
                    continue;

                if(exon.Start > codingStop)
                    break;

                codingBaseCount += min(codingStop, exon.End) - max(refPosition, exon.Start) + 1;
            }
        }
        else
        {
            int codingStop = transData.CodingStart != null ? transData.CodingStart : transData.TransStart;

            if(codingStop >= refPosition)
                return 0;

            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(refPosition < exon.Start)
                    continue;

                if(exon.End < codingStop)
                    break;

                codingBaseCount += min(refPosition, exon.End) - max(codingStop, exon.Start) + 1;
            }
        }

        return codingBaseCount;
    }


    public static int calcStartCodonBases(final NeoEpitope neData)
    {
        // distance from last upstream base to upstream coding start
        final TranscriptData transData = neData.TransData[FS_UP];

        if(transData.CodingStart == null)
            return -1;

        final List<ExonData> exonDataList = transData.exons();
        int refPosition = neData.position(FS_UP);

        int codingBaseCount = 0;

        if(neData.orientation(FS_UP) == NEG_ORIENT)
        {
            int codingPosition = transData.CodingEnd;

            if(codingPosition < refPosition)
                return -1;

            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if (refPosition > exon.End)
                    continue;

                if(exon.Start > codingPosition)
                    break;

                codingBaseCount += min(codingPosition, exon.End) - max(refPosition, exon.Start) + 1;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                int codingPosition = transData.CodingStart;

                if(codingPosition > refPosition)
                    return -1;

                final ExonData exon = exonDataList.get(i);

                if(refPosition < exon.Start)
                    continue;

                if(exon.End < codingPosition)
                    break;

                codingBaseCount += min(refPosition, exon.End) - max(codingPosition, exon.Start) + 1;
            }
        }

        return codingBaseCount;
    }

    public static int findSkippedExonBoundaries(
            final List<TranscriptData> transDataList, final int[] positionBounds, boolean findExonStart, boolean isAcceptor)
    {
        if(positionBounds[SE_START] >= positionBounds[SE_END])
            return 0;

        if(transDataList.isEmpty())
            return 0;

        final Set<Integer> skippedSites = Sets.newHashSet();

        for(TranscriptData transData : transDataList)
        {
            for(ExonData exon : transData.exons())
            {
                if(isAcceptor && exon.Rank == 1)
                    continue;

                if(findExonStart)
                {
                    if(exon.Start <= positionBounds[SE_START])
                        continue;

                    if(exon.Start >= positionBounds[SE_END])
                        break;

                    skippedSites.add(exon.Start);
                }
                else
                {
                    if(exon.End <= positionBounds[SE_START])
                        continue;

                    if(exon.End >= positionBounds[SE_END])
                        break;

                    skippedSites.add(exon.End);
                }
            }
        }

        return skippedSites.size();
    }

    public static String checkTrimBases(final String bases)
    {
        if(bases.length() < 50)
            return bases;

        return bases.substring(0, 50) + "...";
    }

    public static String getAminoAcids(final String baseString, boolean checkStopCodon)
    {
        if(baseString.length() < 3)
            return "";

        String aminoAcidStr = "";
        int index = 0;
        while(index <= baseString.length() - 3)
        {
            String codonBases = baseString.substring(index, index + 3);

            if(checkStopCodon && isStopCodon(codonBases))
            {
                aminoAcidStr += STOP_SYMBOL;
                break;
            }

            String aminoAcid = convertDnaCodonToAminoAcid(codonBases);

            aminoAcidStr += aminoAcid;
            index += 3;
        }

        return aminoAcidStr;
    }

    public static List<PeptideData> generatePeptides(
            final String upAAs, final String novelAAs, final String downAAs, final int[] peptideLengths, int flankingBases)
    {
        final Set<String> uniquePeptides = Sets.newHashSet();
        final List<PeptideData> peptides = Lists.newArrayList();

        // replace any mis-classified upstream stop-codons with 'U' = selenocysteine
        final String fullAminoAcids = upAAs.replaceAll(STOP_SYMBOL, AA_SELENOCYSTEINE) + novelAAs + downAAs;
        int upstreamLength = upAAs.length();
        int novelLength = novelAAs.length();
        boolean hasNovelAAs = !novelAAs.isEmpty();

        int fullLength = fullAminoAcids.length();

        if(fullAminoAcids.endsWith(STOP_SYMBOL))
            --fullLength;

        // eg upstream length = 10, index 0 -> 9, novel at index 10 and downstream at index 11, or no novel and downstream at index 10
        for(int length = peptideLengths[SE_START]; length <= peptideLengths[SE_END]; ++length)
        {
            int currentIndex = 0;

            while(true)
            {
                // skip to point where novel or downstream is included
                if(currentIndex + length <= upstreamLength)
                {
                    ++currentIndex;
                    continue;
                }

                // break if into the downstream (canonical) sequence
                if((!hasNovelAAs && currentIndex >= upstreamLength) || currentIndex >= upstreamLength + novelLength)
                    break;

                if(currentIndex + length > fullLength)
                    break;

                String peptide = fullAminoAcids.substring(currentIndex, currentIndex + length);

                if(!uniquePeptides.contains(peptide))
                {
                    uniquePeptides.add(peptide);

                    String upFlank = "";
                    String downFlank = "";

                    if(flankingBases > 0)
                    {
                        if(currentIndex > 0)
                            upFlank = fullAminoAcids.substring(max(currentIndex - flankingBases, 0), currentIndex);

                        if(currentIndex + length < fullAminoAcids.length())
                        {
                            downFlank = fullAminoAcids.substring(
                                    currentIndex + length, min(currentIndex + length + flankingBases, fullAminoAcids.length()));

                            downFlank = downFlank.replaceAll(STOP_SYMBOL, "");
                        }
                    }

                    peptides.add(new PeptideData(peptide, upFlank, downFlank));
                }

                ++currentIndex;
            }
        }

        return peptides;
    }

    public static String convertHlaTypeForPredictions(final String hlaType)
    {
        // accepted: A*:01:01:01:01 or A*:11:01 or A*11:01 (ie missing the first ':'
        // output: HLA-A0101
        if(hlaType.startsWith("HLA-") && hlaType.length() == 9)
            return hlaType;

        if(hlaType.contains("*") && hlaType.contains(":"))
        {
            String[] geneTypes = hlaType.split("\\*");

            if(geneTypes.length != 2)
                return null;

            String[] types = geneTypes[1].split(":");

            if(types.length == 2 && !types[0].isEmpty() && !types[1].isEmpty()) // like A*11:01
                return "HLA-" + geneTypes[0] + types[0] + types[1];

            if(types.length <= 2 || types[1].isEmpty() || types[2].isEmpty())
                return null;

            return "HLA-" + geneTypes[0] + types[1] + types[2]; // like A*:11:01
        }

        return null;
    }

}
