package com.hartwig.hmftools.linx.gene;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.gene.BreakendTransData.POST_CODING_PHASE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvVarData;

public class BreakendGenePrep
{
    public static void setSvGeneData(
            final List<SvVarData> svList, final EnsemblDataCache ensemblDataCache, boolean applyPromotorDistance, boolean loadBreakendGenes)
    {
        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        if (loadBreakendGenes)
        {
            // only load transcript info for the genes covered
            final List<String> restrictedGeneIds = Lists.newArrayList();

            for (final SvVarData var : svList)
            {
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (be == SE_END && var.isSglBreakend())
                    {
                        // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                        for(final SglMapping mapping : var.getSglMappings())
                        {
                            ensemblDataCache.populateGeneIdList(restrictedGeneIds, mapping.Chromosome, mapping.Position, upstreamDistance);
                        }
                    }
                    else
                    {
                        boolean isStart = isStart(be);
                        ensemblDataCache.populateGeneIdList(restrictedGeneIds, var.chromosome(isStart), var.position(isStart), upstreamDistance);
                    }
                }
            }

            ensemblDataCache.getAlternativeGeneData().stream().filter(x -> !restrictedGeneIds.contains(x.GeneId))
                    .forEach(x -> restrictedGeneIds.add(x.GeneId));

            ensemblDataCache.loadTranscriptData(restrictedGeneIds);
        }

        // associate breakends with transcripts
        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);
                final List<BreakendGeneData> genesList = var.getGenesList(isStart);

                if (be == SE_END && var.isSglBreakend())
                {
                    // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                    for(final SglMapping mapping : var.getSglMappings())
                    {
                        final List<BreakendGeneData> mappingGenes = findGeneAnnotationsBySv(
                                ensemblDataCache, var.id(), isStart, mapping.Chromosome, mapping.Position, mapping.Orientation,
                                upstreamDistance);

                        mappingGenes.forEach(x -> x.setType(var.type()));

                        genesList.addAll(mappingGenes);
                    }
                }
                else
                {
                    genesList.addAll(findGeneAnnotationsBySv(
                            ensemblDataCache, var.id(), isStart, var.chromosome(isStart), var.position(isStart), var.orientation(isStart),
                            upstreamDistance));

                    for (BreakendGeneData gene : genesList)
                    {
                        gene.setSvData(var.getSvData(), var.jcn());
                    }
                }
            }
        }
    }

    public static List<BreakendGeneData> findGeneAnnotationsBySv(
            final EnsemblDataCache ensemblDataCache, int svId, boolean isStart, final String chromosome, int position, byte orientation,
            int upstreamDistance)
    {
        List<BreakendGeneData> geneAnnotations = Lists.newArrayList();

        final List<GeneData> matchedGenes = ensemblDataCache.findGeneRegions(chromosome, position, upstreamDistance);

        // now look up relevant transcript and exon information
        for(final GeneData geneData : matchedGenes)
        {
            final List<TranscriptData> transcriptDataList = ensemblDataCache.getTranscriptDataMap().get(geneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            BreakendGeneData currentGene = new BreakendGeneData(svId, isStart, geneData);

            currentGene.setPositionalData(chromosome, position, orientation);

            // collect up all the relevant exons for each unique transcript to analyse as a collection
            for(TranscriptData transData : transcriptDataList)
            {
                BreakendTransData transcript = createBreakendTranscriptData(transData, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);

                    setAlternativeTranscriptPhasings(transcript, transData.exons(), position, orientation);

                    // annotate with preceding gene info if the up distance isn't set
                    if(!transcript.hasPrevSpliceAcceptorDistance())
                    {
                        setPrecedingGeneDistance(ensemblDataCache, transcript, position);
                    }
                }
            }

            if(currentGene.transcripts().isEmpty() && ensemblDataCache.hasDownstreamGeneAnnotation(geneData))
            {
                // generate a canonical transcript record for the downstream SV position
                final TranscriptData transData = transcriptDataList.stream().filter(x -> x.IsCanonical).findAny().orElse(null);

                if(transData != null)
                {
                    final CodingBaseData cbData = calcCodingBases(transData, position);

                    final BreakendTransData postGeneTrans = new BreakendTransData(
                            currentGene, transData, 1, 1, PHASE_NONE, PHASE_NONE, cbData.CodingBases, cbData.TotalCodingBases);

                    currentGene.addTranscript(postGeneTrans);
                }
            }

            geneAnnotations.add(currentGene);
        }

        final List<GeneData> altMappingGenes = ensemblDataCache.getAlternativeGeneData().stream()
                .filter(x -> x.Chromosome.equals(chromosome))
                .filter(x -> positionWithin(position, x.GeneStart, x.GeneEnd))
                .collect(Collectors.toList());

        for(final GeneData altGeneData : altMappingGenes)
        {
            final List<TranscriptData> transcriptDataList = ensemblDataCache.getTranscriptDataMap().get(altGeneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            final TranscriptData trans = transcriptDataList.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

            BreakendGeneData geneAnnotation = new BreakendGeneData(svId, isStart, altGeneData);

            byte downstreamOrient = altGeneData.Strand == POS_ORIENT ? NEG_ORIENT : POS_ORIENT;
            geneAnnotation.setPositionalData(chromosome, position, downstreamOrient);

            if(trans != null)
            {
                final CodingBaseData cbData = calcCodingBases(trans, position);

                final BreakendTransData altGeneTrans = new BreakendTransData(
                        geneAnnotation, trans, 1, 1, PHASE_NONE, PHASE_NONE, cbData.CodingBases, cbData.TotalCodingBases);

                geneAnnotation.addTranscript(altGeneTrans);
            }

            geneAnnotations.add(geneAnnotation);
        }

        return geneAnnotations;
    }

    public static void setAlternativeTranscriptPhasings(
            final BreakendTransData transcript, final List<ExonData> exonDataList, int position, byte orientation)
    {
        // collect exon phasings before the position on the upstream and after it on the downstream
        boolean isUpstream = (transcript.gene().strand() * orientation) > 0;
        boolean forwardStrand = transcript.TransData.posStrand();

        Map<Integer,Integer> alternativePhasing = transcript.getAlternativePhasing();
        alternativePhasing.clear();

        int transPhase = isUpstream ? transcript.Phase : transcript.Phase;
        int transRank = isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream;

        for (ExonData exonData : exonDataList)
        {
            if(isUpstream == forwardStrand)
            {
                if(exonData.Start > position || transRank == exonData.Rank)
                    break;
            }
            else
            {
                if(position > exonData.End || transRank == exonData.Rank)
                    continue;
            }

            int exonPhase = isUpstream ? exonData.PhaseEnd : exonData.PhaseStart;

            if(exonPhase == PHASE_NONE && !transcript.TransData.nonCoding())
            {
                // switch to -2 as per the 3'UTR convention
                if(forwardStrand && exonData.End > transcript.TransData.CodingEnd)
                    exonPhase = POST_CODING_PHASE;
                else if(!forwardStrand && exonData.Start < transcript.TransData.CodingStart)
                    exonPhase = POST_CODING_PHASE;
            }

            int exonsSkipped;

            if(isUpstream)
            {
                exonsSkipped = max(transRank - exonData.Rank, 0);
            }
            else
            {
                exonsSkipped = max(exonData.Rank - transRank, 0);
            }

            if(exonPhase != transPhase)
            {
                if(isUpstream == forwardStrand)
                {
                    // take the closest to the position
                    alternativePhasing.put(exonPhase, exonsSkipped);
                }
                else
                {
                    // take the first found
                    if(!alternativePhasing.containsKey(exonPhase))
                        alternativePhasing.put(exonPhase, exonsSkipped);
                }
            }
        }
    }

    public static BreakendTransData createBreakendTranscriptData(
            final TranscriptData transData, int position, final BreakendGeneData geneAnnotation)
    {
        final List<ExonData> exonList = transData.exons();

        if(exonList.isEmpty())
            return null;

        boolean isForwardStrand = geneAnnotation.strand() == POS_STRAND;
        boolean isUpstream = geneAnnotation.isUpstream();

        int upExonRank = -1;
        int downExonRank = -1;
        int nextUpDistance = -1;
        int nextDownDistance = -1;
        boolean isCodingTypeOverride = false;
        int phase = PHASE_NONE;

        // first check for a position outside the exon boundaries
        final ExonData firstExon = exonList.get(0);
        final ExonData lastExon = exonList.get(exonList.size()-1);

        // for forward-strand transcripts the current exon is downstream, the previous is upstream
        // and the end-phase is taken from the upstream previous exon, the phase from the current downstream exon

        // for reverse-strand transcripts the current exon is upstream, the previous is downstream
        // and the end-phase is taken from the upstream (current) exon, the phase from the downstream (previous) exon

        // for each exon, the 'phase' is always the phase at the start of the exon in the direction of transcription
        // regardless of strand direction, and 'end_phase' is the phase at the end of the exon

        if(position < firstExon.Start)
        {
            if(isForwardStrand)
            {
                // proceed to the next exon assuming its splice acceptor is required
                final ExonData firstSpaExon = exonList.size() > 1 ? exonList.get(1) : firstExon;
                downExonRank = firstSpaExon.Rank;
                nextDownDistance = firstSpaExon.Start - position;

                isCodingTypeOverride = transData.CodingStart != null && firstSpaExon.Start > transData.CodingStart;

                if(transData.CodingStart != null)
                {
                    if(firstSpaExon.Start > transData.CodingStart)
                        isCodingTypeOverride = true;

                    if(firstSpaExon.Start == transData.CodingStart)
                        phase = PHASE_NONE;
                    else
                        phase = firstSpaExon.PhaseStart;
                }

                upExonRank = 0;
            }
            else
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position post-coding
                return null;
            }
        }
        else if(position > lastExon.End)
        {
            if(!isForwardStrand)
            {
                final ExonData firstSpaExon = exonList.size() > 1 ? exonList.get(exonList.size()-2) : lastExon;
                downExonRank = firstSpaExon.Rank;
                nextDownDistance = position - lastExon.End;

                if(transData.CodingEnd != null)
                {
                    if(firstSpaExon.End < transData.CodingEnd)
                        isCodingTypeOverride = true;

                    if(firstSpaExon.End == transData.CodingEnd)
                        phase = PHASE_NONE;
                    else
                        phase = firstSpaExon.PhaseStart;
                }

                upExonRank = 0;
            }
            else
            {
                return null;
            }
        }
        else
        {
            for (int index = 0; index < exonList.size(); ++index)
            {
                final ExonData exonData = exonList.get(index);

                if (positionWithin(position, exonData.Start, exonData.End))
                {
                    // falls within an exon
                    upExonRank = downExonRank = exonData.Rank;

                    // set distance to next and previous splice acceptor
                    if(isForwardStrand)
                    {
                        nextUpDistance = position - exonData.Start;

                        if(index < exonList.size() - 1)
                        {
                            final ExonData nextExonData = exonList.get(index + 1);
                            nextDownDistance = nextExonData.Start - position;
                        }
                    }
                    else
                    {
                        nextUpDistance = exonData.End - position;

                        if(index > 1)
                        {
                            // first splice acceptor is the second exon (or later on)
                            final ExonData prevExonData = exonList.get(index - 1);
                            nextDownDistance = position - prevExonData.End;
                        }
                    }

                    phase = isUpstream ? exonData.PhaseStart : exonData.PhaseEnd;

                    break;
                }
                else if(position < exonData.Start)
                {
                    // position falls between this exon and the previous one
                    final ExonData prevExonData = exonList.get(index-1);

                    if(isForwardStrand)
                    {
                        // the current exon is downstream, the previous one is upstream
                        upExonRank = prevExonData.Rank;
                        downExonRank = exonData.Rank;
                        nextDownDistance = exonData.Start - position;
                        nextUpDistance = position - prevExonData.End;
                    }
                    else
                    {
                        // the current exon is earlier in rank
                        // the previous exon in the list has the higher rank and is downstream
                        // the start of the next exon (ie previous here) uses 'phase' for the downstream as normal
                        upExonRank = exonData.Rank;
                        downExonRank = prevExonData.Rank;
                        nextUpDistance = exonData.Start - position;
                        nextDownDistance = position - prevExonData.End;
                    }

                    if(isUpstream)
                    {
                        if(isForwardStrand)
                            phase = prevExonData.PhaseEnd;
                        else
                            phase = exonData.PhaseEnd;
                    }
                    else
                    {
                        // if coding starts on the first base of the next exon, use -1
                        if(isForwardStrand)
                        {
                            if(transData.CodingStart != null && transData.CodingStart == exonData.Start)
                                phase = PHASE_NONE;
                            else
                                phase = exonData.PhaseStart;
                        }
                        else
                        {
                            if(transData.CodingEnd != null && transData.CodingEnd == prevExonData.End)
                                phase = PHASE_NONE;
                            else
                                phase = prevExonData.PhaseStart;
                        }
                    }

                    break;
                }
            }
        }

        // now calculate coding bases for this transcript
        // for the given position, determine how many coding bases occur prior to the position
        // in the direction of the transcript

        final CodingBaseData cbData = calcCodingBases(transData, position);

        BreakendTransData transcript = new BreakendTransData(geneAnnotation, transData,
                upExonRank, downExonRank, phase, cbData.Phase, cbData.CodingBases, cbData.TotalCodingBases);

        // if not set, leave the previous exon null and it will be taken from the closest upstream gene
        transcript.setSpliceAcceptorDistance(true, nextUpDistance >= 0 ? nextUpDistance : null);
        transcript.setSpliceAcceptorDistance(false, nextDownDistance >= 0 ? nextDownDistance : null);

        if(isCodingTypeOverride)
            transcript.setCodingType(CODING);

        return transcript;
    }

    private static void setPrecedingGeneDistance(final EnsemblDataCache ensemblDataCache, final BreakendTransData transcript, int position)
    {
        // annotate with preceding gene info if the up distance isn't set
        int precedingGeneSAPos = ensemblDataCache.findPrecedingGeneSpliceAcceptorPosition(transcript.transId());

        if(precedingGeneSAPos >= 0)
        {
            // if the breakend is after (higher for +ve strand) the nearest preceding splice acceptor, then the distance will be positive
            // and mean that the transcript isn't interrupted when used in a downstream fusion
            int preDistance = transcript.gene().strand() == POS_STRAND ? position - precedingGeneSAPos : precedingGeneSAPos - position;
            transcript.setSpliceAcceptorDistance(true, preDistance);
        }
    }

}
