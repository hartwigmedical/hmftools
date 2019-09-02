package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.linx.types.SvConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter.GENE_TYPE_PSEUDOGENE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class PseudoGeneFinder
{
    private final VisualiserWriter mVisWriter;
    private SvGeneTranscriptCollection mGeneTransCache;

    private static final Logger LOGGER = LogManager.getLogger(PseudoGeneFinder.class);

    public PseudoGeneFinder(final VisualiserWriter visWriter)
    {
        mVisWriter = visWriter;
        mGeneTransCache = null;
    }

    public void setGeneTransCache(final SvGeneTranscriptCollection geneTransCache)
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

                        Map<Integer,int[]> exonPositionOffsets = Maps.newHashMap();

                        final int selectedTransId = maxTrans.TransId;

                        // cache this info against the link so it can be recorded in output files and the DB
                        for(final SvLinkedPair pair : matchedPairs)
                        {
                            final PseudoGeneMatch pseudoMatch = pairMatchesMap.get(pair).stream()
                                    .filter(x -> x.TransId == selectedTransId).findFirst().orElse(null);

                            if(pseudoMatch != null)
                            {
                                int[] exonPosOffsets = new int[SE_PAIR];

                                String exonMatchData = String.format("%s;%s;%d;%d",
                                        gene.GeneName, maxTrans.TransName, pseudoMatch.ExonRank, pseudoMatch.ExonLength);

                                for(int se = SE_START; se <= SE_END; ++se)
                                {
                                    final SvBreakend breakend = pair.getBreakend(se);

                                    if (pseudoMatch.isHomologyMatch(isStart(se)))
                                    {
                                        exonPosOffsets[se] = pseudoMatch.HomologyOffset[se];
                                    }

                                    boolean homMismatch = hasHomologyMismatch(breakend.getSV(), pairMatchesMap, selectedTransId);

                                    exonMatchData += String.format(";%d;%d;%s",
                                            pseudoMatch.HomologyOffset[se], pseudoMatch.PositionMismatch[se], homMismatch);
                                }

                                pair.setExonMatchData(exonMatchData);

                                exonPositionOffsets.put(pseudoMatch.ExonRank, exonPosOffsets);
                            }
                        }

                        mVisWriter.addGeneExonData(cluster.id(), gene.StableId, gene.GeneName,
                                maxTrans.TransName, maxTrans.TransId, gene.chromosome(), GENE_TYPE_PSEUDOGENE, exonPositionOffsets);
                    }
                }
            }
        }
    }

    private List<PseudoGeneMatch> findPseudoGeneExonMatches(
            final GeneAnnotation gene, long posStart, long posEnd, int startHomologyLength, int endHomologyLength)
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

        return homOffsetStart == homOffsetEnd;
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
}
