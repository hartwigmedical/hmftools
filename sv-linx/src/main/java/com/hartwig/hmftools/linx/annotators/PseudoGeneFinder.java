package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.linx.types.SvConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter.GENE_TYPE_PSEUDOGENE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneData;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;

public class PseudoGeneFinder
{
    private final VisualiserWriter mVisWriter;
    private EnsemblDataCache mGeneTransCache;

    public PseudoGeneFinder(final VisualiserWriter visWriter)
    {
        mVisWriter = visWriter;
        mGeneTransCache = null;
    }

    public void setGeneTransCache(final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;
    }

    public void checkPseudoGeneAnnotations(final List<SvCluster> clusters)
    {
        if(mGeneTransCache == null)
            return;

        for(final SvCluster cluster : clusters)
        {
            Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap = Maps.newHashMap();
            List<GeneAnnotation> matchedGenes = Lists.newArrayList();

            for(final SvLinkedPair pair : cluster.getLinkedPairs())
            {
                if(pair.length() > SHORT_TI_LENGTH * 8)
                    continue;

                final SvBreakend lower = pair.getBreakend(true);
                final SvBreakend upper = pair.getBreakend(false);

                // for any TI falling within the same gene, check for an exon boundary match
                if(lower.getSV().getGenesList(lower.usesStart()).isEmpty() || upper.getSV().getGenesList(upper.usesStart()).isEmpty())
                    continue;

                final String lowerHomology = lower.usesStart() ?
                        lower.getSV().getSvData().startHomologySequence() : lower.getSV().getSvData().endHomologySequence();

                final String upperHomology = lower.usesStart() ?
                        upper.getSV().getSvData().startHomologySequence() : upper.getSV().getSvData().endHomologySequence();

                for(final GeneAnnotation gene1 : lower.getSV().getGenesList(lower.usesStart()))
                {
                    for(final GeneAnnotation gene2 : upper.getSV().getGenesList(upper.usesStart()))
                    {
                        if(!gene1.GeneName.equals(gene2.GeneName))
                            continue;

                        final List<PseudoGeneMatch> pseudoMatches = findPseudoGeneExonMatches(
                                gene1, lower.position(), upper.position(), lowerHomology.length(), upperHomology.length());

                        if(!pseudoMatches.isEmpty())
                        {
                            pairMatchesMap.put(pair, pseudoMatches);

                            if(!matchedGenes.stream().anyMatch(x -> x.GeneName.matches(gene1.GeneName)))
                            {
                                matchedGenes.add(gene1);
                            }
                        }
                    }
                }
            }

            if(!pairMatchesMap.isEmpty())
            {
                // select the most common transcript and report it for visualisation
                for (final GeneAnnotation gene : matchedGenes)
                {
                    // find the most frequent transcript
                    final PseudoGeneMatch maxTrans = findMostCommonTranscript(gene.GeneName, pairMatchesMap);

                    if(maxTrans != null)
                    {
                        // mark a pseudogene link match if either both breakends are an exact match
                        // or one of the SVs also matches another pseudogene
                        List<SvLinkedPair> matchedPairs = pairMatchesMap.keySet().stream()
                                .filter(x -> isRelevantMatch(x, pairMatchesMap, maxTrans.TransId))
                                .collect(Collectors.toList());

                        if(matchedPairs.isEmpty())
                            continue;

                        VisGeneData geneData = new VisGeneData(cluster.id(), gene.StableId, gene.GeneName,
                                maxTrans.TransName, maxTrans.TransId, gene.chromosome(), GENE_TYPE_PSEUDOGENE);

                        final int selectedTransId = maxTrans.TransId;

                        // cache this info against the link so it can be recorded in output files and the DB
                        for(final SvLinkedPair pair : matchedPairs)
                        {
                            final PseudoGeneMatch pseudoMatch = pairMatchesMap.get(pair).stream()
                                    .filter(x -> x.TransId == selectedTransId).findFirst().orElse(null);

                            if(pseudoMatch != null)
                            {
                                int[] exonPosOffsets = new int[SE_PAIR];
                                int[] exonsLost = new int[SE_PAIR];

                                String exonMatchData = String.format("%s;%s;%d;%d",
                                        gene.GeneName, maxTrans.TransName, pseudoMatch.ExonRank, pseudoMatch.ExonLength);

                                for(int se = SE_START; se <= SE_END; ++se)
                                {
                                    final SvBreakend breakend = pair.getBreakend(se);

                                    if (pseudoMatch.isHomologyMatch(isStart(se)))
                                    {
                                        exonPosOffsets[se] = pseudoMatch.HomologyOffset[se];
                                    }
                                    else if(pseudoMatch.PositionMismatch[se] < 0)
                                    {
                                        exonPosOffsets[se] = 0;
                                        exonsLost[se] = pseudoMatch.PositionMismatch[se];
                                    }

                                    boolean homMismatch = hasHomologyMismatch(breakend.getSV(), pairMatchesMap, selectedTransId);

                                    exonMatchData += String.format(";%d;%d;%s",
                                            pseudoMatch.HomologyOffset[se], pseudoMatch.PositionMismatch[se], homMismatch);
                                }

                                pair.setExonMatchData(exonMatchData);

                                geneData.ExonPositionOffsets.put(pseudoMatch.ExonRank, exonPosOffsets);

                                if(exonsLost[SE_START] != 0 || exonsLost[SE_END] != 0)
                                    geneData.ExonsLostOffsets.put(pseudoMatch.ExonRank, exonsLost);
                            }
                        }

                        final TranscriptData transData = mGeneTransCache.getTranscriptData(geneData.GeneId, geneData.TransName);

                        if(transData == null || transData.exons().isEmpty())
                            continue;

                        for (final ExonData exonData : transData.exons())
                        {
                            boolean hasPosOffsets = geneData.ExonPositionOffsets.containsKey(exonData.ExonRank);

                            int[] exonsLost = geneData.ExonsLostOffsets.get(exonData.ExonRank);

                            if(hasPosOffsets && exonsLost == null)
                                continue;

                            // every exon will have an entry to show the full set of exons
                            if(!hasPosOffsets)
                            {
                                geneData.ExonPositionOffsets.put(exonData.ExonRank, new int[] { 0, 0 });
                            }

                            // now work out position adjustments for the exons lost
                            if(exonsLost == null)
                            {
                                geneData.ExonsLostOffsets.put(exonData.ExonRank, new int[] { 0, 0});
                            }
                            else
                            {
                                int exonLength = exonData.ExonEnd - exonData.ExonStart;

                                // if say X bases have been lost from the start, then set the end to factor this in
                                if(exonsLost[SE_START] < 0)
                                {
                                    exonsLost[SE_END] = -(int)(exonLength + exonsLost[SE_START]);
                                    exonsLost[SE_START] = 0;
                                }
                                else if(exonsLost[SE_END] < 0)
                                {
                                    exonsLost[SE_START] = (int)(exonLength + exonsLost[SE_END]);
                                    exonsLost[SE_END] = 0;
                                }
                            }
                        }

                        mVisWriter.addGeneExonData(geneData);
                    }
                }
            }
        }
    }

    private List<PseudoGeneMatch> findPseudoGeneExonMatches(
            final GeneAnnotation gene, int posStart, int posEnd, int startHomologyLength, int endHomologyLength)
    {
        List<PseudoGeneMatch> pseudoMatches = Lists.newArrayList();

        List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(gene.StableId);

        for(final TranscriptData transData : transDataList)
        {
            int exonCount = transData.exons().size();

            for (int i = 0; i < exonCount; ++i)
            {
                final ExonData exonData = transData.exons().get(i);

                boolean startWithinHomology = abs(exonData.ExonStart - posStart) <= startHomologyLength;
                boolean endWithinHomology = abs(exonData.ExonEnd - posEnd) <= endHomologyLength;

                if(!startWithinHomology && !endWithinHomology)
                    continue;

                // skip if the non-matching end is outside the transcript
                if(!startWithinHomology && posStart < transData.TransStart)
                    continue;

                if(!endWithinHomology && posEnd > transData.TransEnd)
                    continue;

                PseudoGeneMatch pseudoMatch = new PseudoGeneMatch(
                        gene.GeneName, transData.TransId, transData.TransName, exonData.ExonRank, exonData.length());

                if(startWithinHomology)
                {
                    pseudoMatch.HomologyOffset[SE_START] = (int)(posStart - exonData.ExonStart);
                }
                else
                {
                    // record a position within the exon as negative
                    pseudoMatch.PositionMismatch[SE_START] = (int)(exonData.ExonStart - posStart);
                }

                if(endWithinHomology)
                {
                    pseudoMatch.HomologyOffset[SE_END] = (int)(posEnd - exonData.ExonEnd);
                }
                else
                {
                    pseudoMatch.PositionMismatch[SE_END] = (int)(posEnd - exonData.ExonEnd);
                }

                pseudoMatches.add(pseudoMatch);
            }
        }

        return pseudoMatches;
    }

    private boolean hasHomologyMismatch(final SvVarData var, final Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap, int selectedTransId)
    {
        if(var.isSglBreakend())
            return false;

        final Integer homOffsetStart = getHomologyOffset(var.getBreakend(true), pairMatchesMap, selectedTransId);
        final Integer homOffsetEnd = getHomologyOffset(var.getBreakend(false), pairMatchesMap, selectedTransId);

        if(homOffsetStart == null || homOffsetEnd == null)
            return false;

        return homOffsetStart != homOffsetEnd;
    }

    private Integer getHomologyOffset(
            final SvBreakend breakend, final Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap, int selectedTransId)
    {
        for (Map.Entry<SvLinkedPair, List<PseudoGeneMatch>> entry : pairMatchesMap.entrySet())
        {
            final SvLinkedPair pair = entry.getKey();

            if (pair.hasBreakend(breakend))
            {
                final PseudoGeneMatch pseudoMatch = pairMatchesMap.get(pair).stream()
                        .filter(x -> x.TransId == selectedTransId).findFirst().orElse(null);

                if(pseudoMatch != null)
                {
                    if(pair.getBreakend(true) == breakend)
                        return pseudoMatch.isHomologyMatch(true) ? pseudoMatch.HomologyOffset[SE_START] : null;
                    else
                        return pseudoMatch.isHomologyMatch(false) ? pseudoMatch.HomologyOffset[SE_END] : null;
                }
            }
        }

        return null;
    }

    private PseudoGeneMatch findPseudoMatch(
            final SvBreakend breakend, final Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap, int selectedTransId)
    {
        for (Map.Entry<SvLinkedPair, List<PseudoGeneMatch>> entry : pairMatchesMap.entrySet())
        {
            final SvLinkedPair pair = entry.getKey();

            if (pair.hasBreakend(breakend))
            {
                final PseudoGeneMatch pseudoMatch = pairMatchesMap.get(pair).stream()
                        .filter(x -> x.TransId == selectedTransId).findFirst().orElse(null);

                return pseudoMatch;
            }
        }

        return null;
    }

    private boolean isRelevantMatch(final SvLinkedPair pair, final Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap, int selectedTransId)
    {
        final PseudoGeneMatch pseudoMatch = pairMatchesMap.get(pair).stream()
                .filter(x -> x.TransId == selectedTransId).findFirst().orElse(null);

        if(pseudoMatch == null)
            return false;

        if(pseudoMatch.isHomologyMatch(true) && pseudoMatch.isHomologyMatch(false))
            return true;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final SvBreakend breakend = pair.getBreakend(se);
            final SvBreakend otherBreakend = breakend.getOtherBreakend();

            if(otherBreakend == null)
                continue;

            if(findPseudoMatch(otherBreakend, pairMatchesMap, selectedTransId) != null)
                return true;
        }

        return false;
    }

    private PseudoGeneMatch findMostCommonTranscript(final String geneName, final Map<SvLinkedPair,List<PseudoGeneMatch>> pairMatchesMap)
    {
        Map<Integer, Integer> transcriptMatches = Maps.newHashMap();
        PseudoGeneMatch maxTrans = null;
        int maxTransIdCount = 0;

        for (Map.Entry<SvLinkedPair, List<PseudoGeneMatch>> entry : pairMatchesMap.entrySet())
        {
            for (PseudoGeneMatch pseudoMatch : entry.getValue())
            {
                if(!pseudoMatch.Gene.equals(geneName))
                    continue;

                int count = transcriptMatches.containsKey(pseudoMatch.TransId) ? transcriptMatches.get(pseudoMatch.TransId) : 0;
                ++count;

                transcriptMatches.put(pseudoMatch.TransId, count);

                if(count > maxTransIdCount)
                {
                    maxTransIdCount = count;
                    maxTrans = pseudoMatch;
                }
            }
        }

        return maxTrans;
    }

    public boolean variantMatchesPseudogeneExons(final SvVarData var)
    {
        if(mGeneTransCache == null)
            return false;

        if(var.isSglBreakend())
            return false;

        final List<GeneAnnotation> genesStart = var.getGenesList(true);
        final List<GeneAnnotation> genesEnd = var.getGenesList(true);

        if(genesStart.isEmpty() || genesEnd.isEmpty())
            return false;

        int posStart = var.position(true);
        int posEnd = var.position(false);

        int startHomologyLength = var.getSvData().startHomologySequence().length();
        int endHomologyLength = var.getSvData().endHomologySequence().length();

        for(final GeneAnnotation gene : genesStart)
        {
            if(!genesEnd.stream().anyMatch(x -> x.GeneName.equals(gene.GeneName)))
                continue;

            List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(gene.StableId);

            for(final TranscriptData transData : transDataList)
            {
                if(!transData.exons().stream().anyMatch(x ->
                        abs(x.ExonStart - posStart) <= startHomologyLength || abs(x.ExonEnd - posStart) <= startHomologyLength))
                {
                    continue;
                }

                if(transData.exons().stream().anyMatch(x ->
                        abs(x.ExonStart - posEnd) <= endHomologyLength || abs(x.ExonEnd - posEnd) <= endHomologyLength))
                {
                    return true;
                }
            }
        }

        return false;
    }

}
