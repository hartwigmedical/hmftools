package com.hartwig.hmftools.svannotation;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_5P_UTR;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_0;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_1;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_CODING_2;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.GENE_PHASING_REGION_MAX;
import static com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData.phaseToRegion;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.TRANS_CODING_TYPE_CODING;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.validFusionTranscript;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<String, List<TranscriptExonData>> mGeneTransExonDataMap;// keyed by GeneId (aka StableId)
    private Map<String, List<EnsemblGeneData>> mChromosomeGeneDataMap;
    private Map<String, List<EnsemblGeneData>> mChromosomeReverseGeneDataMap; // order by gene end not start
    private Map<Integer, List<TranscriptProteinData>> mEnsemblProteinDataMap;

    private BufferedWriter mBreakendWriter;

    // to get a wider range of candidate genes and filter by promotor distance later on
    public static int PRE_GENE_PROMOTOR_DISTANCE = 100000;

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mGeneTransExonDataMap = Maps.newHashMap();
        mChromosomeGeneDataMap = Maps.newHashMap();
        mChromosomeReverseGeneDataMap = Maps.newHashMap();
        mEnsemblProteinDataMap = Maps.newHashMap();
        mBreakendWriter = null;
    }

    public void setDataPath(final String dataPath)
    {
        mDataPath = dataPath;

        if(!mDataPath.endsWith(File.separator))
            mDataPath += File.separator;
    }

    public boolean hasCachedEnsemblData()
    {
        return !mChromosomeGeneDataMap.isEmpty() && !mGeneTransExonDataMap.isEmpty();
    }

    public final Map<String, List<TranscriptExonData>> getGeneExonDataMap() { return mGeneTransExonDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrGeneDataMap() { return mChromosomeGeneDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrReverseGeneDataMap() { return mChromosomeReverseGeneDataMap; }
    public Map<Integer, List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }

    public final EnsemblGeneData getGeneDataByName(final String geneName)
    {
        return getGeneData(geneName, true);
    }

    public final EnsemblGeneData getGeneDataById(final String geneId)
    {
        return getGeneData(geneId, false);
    }

    private final EnsemblGeneData getGeneData(final String gene, boolean byName)
    {
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChromosomeGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                if((byName && geneData.GeneName.equals(gene)) || (!byName && geneData.GeneId.equals(gene)))
                    return geneData;
            }
        }

        return null;
    }

    public List<TranscriptExonData> getTransExonData(final String geneId)
    {
        return mGeneTransExonDataMap.get(geneId);
    }

    public void populateGeneIdList(Map<String,Boolean> uniqueGeneIds, final String chromosome, long position, int upstreamDistance)
    {
        // find the unique set of geneIds
        final List<EnsemblGeneData> geneRegions = mChromosomeGeneDataMap.get(chromosome);

        if (geneRegions == null)
            return;

        List<EnsemblGeneData> matchedGenes = findGeneRegions(position, geneRegions, upstreamDistance);

        for (final EnsemblGeneData geneData : matchedGenes)
        {
            uniqueGeneIds.put(geneData.GeneId,true);
        }
    }

    public List<GeneAnnotation> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, long position,
            int upstreamDistance)
    {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        final List<EnsemblGeneData> geneRegions = mChromosomeGeneDataMap.get(chromosome);
        final List<EnsemblGeneData> geneRegionsReversed = mChromosomeReverseGeneDataMap.get(chromosome);

        if(geneRegions == null)
            return geneAnnotations;

        final List<EnsemblGeneData> matchedGenes = findGeneRegions(position, geneRegions, upstreamDistance);

        // now look up relevant transcript and exon information
        for(final EnsemblGeneData geneData : matchedGenes)
        {
            final List<TranscriptExonData> transExonDataList = mGeneTransExonDataMap.get(geneData.GeneId);

            if (transExonDataList == null || transExonDataList.isEmpty())
                continue;

            GeneAnnotation currentGene = new GeneAnnotation(svId, isStart, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.Synonyms, geneData.EntrezIds, geneData.KaryotypeBand);

            currentGene.setGeneData(geneData);

            // collect up all the relevant exons for each unique transcript to analyse as a collection

            int teIndex = 0;
            List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

            while(!transcriptExons.isEmpty())
            {
                Transcript transcript = extractTranscriptExonData(transcriptExons, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);

                    // annotate with preceding gene info if the up distance isn't set
                    if(transcript.exonDistanceUp() == -1)
                    {
                        setPrecedingGeneDistance(transcript, geneData, position, geneRegions, geneRegionsReversed);
                    }
                }

                teIndex += transcriptExons.size();
                transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    private void setPrecedingGeneDistance(
            Transcript transcript, final EnsemblGeneData geneData, long position,
            final List<EnsemblGeneData> geneRegions, final List<EnsemblGeneData> geneRegionsReversed)
    {
        // annotate with preceding gene info if the up distance isn't set
        long precedingGeneSAPos = findPrecedingGeneSpliceAcceptorPosition(geneData, geneData.Strand == 1 ? geneRegionsReversed : geneRegions);

        if(precedingGeneSAPos >= 0)
        {
            long preDistance = geneData.Strand == 1 ? position - precedingGeneSAPos : precedingGeneSAPos - position;
            transcript.setExonDistances((int)preDistance, transcript.exonDistanceDown());
        }
    }

    public static List<TranscriptExonData> nextTranscriptExons(final List<TranscriptExonData> transExonDataList, int currentIndex)
    {
        List<TranscriptExonData> transcriptExons = Lists.newArrayList();

        if(currentIndex >= transExonDataList.size())
            return transcriptExons;

        int transId = transExonDataList.get(currentIndex).TransId;

        int j = currentIndex;
        for(; j < transExonDataList.size(); ++j)
        {
            if(transExonDataList.get(j).TransId != transId)
                break;

            transcriptExons.add(transExonDataList.get(j));
        }

        return transcriptExons;
    }

    public final List<TranscriptExonData> getTranscriptExons(final String geneId, final String transcriptId)
    {
        List<TranscriptExonData> exonDataList = Lists.newArrayList();

        final List<TranscriptExonData> transExonDataList = mGeneTransExonDataMap.get(geneId);

        if (transExonDataList == null || transExonDataList.isEmpty())
            return exonDataList;

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while(!transcriptExons.isEmpty())
        {
            final TranscriptExonData firstExon = transcriptExons.get(0);

            if(transcriptId.isEmpty() && firstExon.IsCanonical)
                return transcriptExons;
            else if(firstExon.TransName.equals(transcriptId))
                return transcriptExons;

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        return exonDataList;
    }

    public final List<EnsemblGeneData> findGenesByRegion(final String chromosome, long posStart, long posEnd)
    {
        List<EnsemblGeneData> genesList = Lists.newArrayList();

        final List<EnsemblGeneData> geneDataList = mChromosomeGeneDataMap.get(chromosome);

        for(final EnsemblGeneData geneData : geneDataList)
        {
            if(posStart <= geneData.GeneStart && posEnd >= geneData.GeneEnd)
            {
                genesList.add(geneData);
            }
        }

        return genesList;
    }

    public final String getCanonicalTranscriptId(final EnsemblGeneData geneData)
    {
        final List<TranscriptExonData> transExonDataList = mGeneTransExonDataMap.get(geneData.GeneId);

        if (transExonDataList == null || transExonDataList.isEmpty())
            return "";

        for (final TranscriptExonData transData : transExonDataList)
        {
            if (transData.IsCanonical)
                return transData.TransName;
        }

        return "";
    }

    private List<EnsemblGeneData> findGeneRegions(long position, List<EnsemblGeneData> geneDataList, int upstreamDistance)
    {
        List<EnsemblGeneData> matchedGenes = Lists.newArrayList();

        for(final EnsemblGeneData geneData : geneDataList)
        {
            long geneStartRange = geneData.Strand == 1 ? geneData.GeneStart - upstreamDistance : geneData.GeneStart;
            long geneEndRange = geneData.Strand == 1 ? geneData.GeneEnd : geneData.GeneEnd + upstreamDistance;

            if(position >= geneStartRange && position <= geneEndRange)
            {
                matchedGenes.add(geneData);
            }
        }

        return matchedGenes;
    }

    private long findPrecedingGeneSpliceAcceptorPosition(final EnsemblGeneData refGene, List<EnsemblGeneData> geneDataList)
    {
        // find the first upstream non-overlapping gene with a splice acceptor
        if(refGene.Strand == 1)
        {
            for(int i = refGene.getReverseListIndex() - 1; i >= 0; --i)
            {
                final EnsemblGeneData gene = geneDataList.get(i);

                if(gene.Strand != refGene.Strand)
                    continue;

                if(gene.GeneEnd > refGene.GeneStart)
                    continue;

                List<TranscriptExonData> transExonData = getTranscriptExons(gene.GeneId, "");

                if(transExonData.size() <= 1)
                    continue;

                // otherwise taken the start of the last exon
                TranscriptExonData lastExonData = transExonData.get(transExonData.size() - 1);
                return lastExonData.ExonStart;
            }
        }
        else
        {
            for(int i = refGene.getListIndex() + 1; i < geneDataList.size(); ++i)
            {
                final EnsemblGeneData gene = geneDataList.get(i);

                if(gene.Strand != refGene.Strand)
                    continue;

                if(gene.GeneStart < refGene.GeneEnd)
                    continue;

                List<TranscriptExonData> transExonData = getTranscriptExons(gene.GeneId, "");

                if(transExonData.size() <= 1)
                    continue;

                // otherwise taken the start of the last exon
                TranscriptExonData lastExonData = transExonData.get(0);
                return lastExonData.ExonEnd;
            }
        }

        return -1;
    }

    public static Transcript extractTranscriptExonData(final List<TranscriptExonData> transcriptExons, long position, final GeneAnnotation geneAnnotation)
    {
        int exonMax = transcriptExons.size();

        final TranscriptExonData first = transcriptExons.get(0);

        boolean isForwardStrand = geneAnnotation.Strand == 1;

        int upExonRank = -1;
        int upExonPhase = -1;
        int downExonRank = -1;
        int downExonPhase = -1;
        long nextUpDistance = -1;
        long nextDownDistance = -1;
        boolean isCodingTypeOverride = false;

        // first check for a position outside the exon boundaries
        final TranscriptExonData firstExon = transcriptExons.get(0);
        final TranscriptExonData lastExon = transcriptExons.get(transcriptExons.size()-1);

        // for forward-strand transcripts the current exon is downstream, the previous is upstream
        // and the end-phase is taken from the upstream previous exon, the phase from the current downstream exon

        // for reverse-strand transcripts the current exon is upstream, the previous is downstream
        // and the end-phase is taken from the upstream (current) exon, the phase from the downstream (previous) exon

        if(position < firstExon.ExonStart)
        {
            if(isForwardStrand)
            {
                // proceed to the next exon assuming its splice acceptor is required
                final TranscriptExonData firstSpaExon = transcriptExons.size() > 1 ? transcriptExons.get(1) : firstExon;
                downExonRank = firstSpaExon.ExonRank;
                downExonPhase = firstSpaExon.ExonPhase;
                nextDownDistance = firstSpaExon.ExonStart - position;

                // correct the phasing if the next exon starts the coding region
                if(firstSpaExon.CodingStart != null && firstSpaExon.ExonStart == firstSpaExon.CodingStart)
                    downExonPhase = -1;

                isCodingTypeOverride = firstSpaExon.CodingStart != null && firstSpaExon.ExonStart > firstSpaExon.CodingStart;

                upExonRank = 0;
                upExonPhase = -1;
            }
            else
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                return null;
            }
        }
        else if(position > lastExon.ExonEnd)
        {
            if(!isForwardStrand)
            {
                final TranscriptExonData firstSpaExon = transcriptExons.size() > 1 ? transcriptExons.get(transcriptExons.size()-2) : lastExon;
                downExonRank = firstSpaExon.ExonRank;
                downExonPhase = firstSpaExon.ExonPhase;
                nextDownDistance = position - lastExon.ExonEnd;

                if(firstSpaExon.CodingEnd != null && firstSpaExon.ExonEnd == firstSpaExon.CodingEnd)
                    downExonPhase = -1;

                isCodingTypeOverride = firstSpaExon.CodingEnd != null && firstSpaExon.ExonEnd < firstSpaExon.CodingEnd;

                upExonRank = 0;
                upExonPhase = -1;
            }
            else
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                return null;
            }
        }
        else
        {
            for (int index = 0; index < transcriptExons.size(); ++index)
            {
                final TranscriptExonData exonData = transcriptExons.get(index);

                if (position >= exonData.ExonStart && position <= exonData.ExonEnd)
                {
                    // falls within an exon
                    upExonRank = downExonRank = exonData.ExonRank;
                    upExonPhase = exonData.ExonPhase;
                    downExonPhase = exonData.ExonPhaseEnd;
                    nextDownDistance = isForwardStrand ? exonData.ExonEnd - position : position - exonData.ExonStart;
                    nextUpDistance = isForwardStrand ? position - exonData.ExonStart : exonData.ExonEnd - position;
                    break;
                }
                else if(position < exonData.ExonStart)
                {
                    // position falls between this exon and the previous one
                    final TranscriptExonData prevExonData = transcriptExons.get(index-1);

                    if(isForwardStrand)
                    {
                        // the current exon is downstream, the prevous one is upstream
                        upExonRank = prevExonData.ExonRank;
                        upExonPhase = prevExonData.ExonPhaseEnd;
                        downExonRank = exonData.ExonRank;
                        downExonPhase = exonData.ExonPhase;
                        nextDownDistance = exonData.ExonStart - position;
                        nextUpDistance = position - prevExonData.ExonEnd;

                        if(exonData.CodingStart != null && exonData.ExonStart == exonData.CodingStart)
                            downExonPhase = -1;
                    }
                    else
                    {
                        // the previous exon in the list has the higher rank and is dowstream
                        upExonRank = exonData.ExonRank;
                        upExonPhase = exonData.ExonPhaseEnd;
                        downExonRank = prevExonData.ExonRank;
                        downExonPhase = prevExonData.ExonPhase;
                        nextUpDistance = exonData.ExonStart - position;
                        nextDownDistance = position - prevExonData.ExonEnd;

                        if(exonData.CodingEnd != null && prevExonData.ExonEnd == exonData.CodingEnd)
                            downExonPhase = -1;
                    }

                    break;
                }
            }
        }

        // now calculate coding bases for this transcript
        // for the given position, determine how many coding bases occur prior to the position
        // in the direction of the transcript

        boolean isCoding = first.CodingStart != null && first.CodingEnd != null;
        long codingStart = first.CodingStart != null ? first.CodingStart : 0;
        long codingEnd = first.CodingEnd != null ? first.CodingEnd : 0;
        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        long codingBases = 0;
        long totalCodingBases = 0;

        if(isCoding)
        {
            for (TranscriptExonData exonData : transcriptExons)
            {
                long exonStart = exonData.ExonStart;
                long exonEnd = exonData.ExonEnd;

                if (!inCodingRegion)
                {
                    if (exonEnd >= codingStart)
                    {
                        // coding region begins in this exon
                        inCodingRegion = true;

                        totalCodingBases += exonEnd - codingStart + 1;

                        // check whether the position falls in this exon and if so before or after the coding start
                        if (position >= codingStart)
                        {
                            if (position < exonEnd)
                                codingBases += position - codingStart + 1;
                            else
                                codingBases += exonEnd - codingStart + 1;
                        }
                    }
                }
                else if (!codingRegionEnded)
                {
                    if (exonStart > codingEnd)
                    {
                        codingRegionEnded = true;
                    }
                    else if (exonEnd > codingEnd)
                    {
                        // coding region ends in this exon
                        codingRegionEnded = true;

                        totalCodingBases += codingEnd - exonStart + 1;

                        if (position >= exonStart)
                        {
                            if (position < codingEnd)
                                codingBases += position - exonStart + 1;
                            else
                                codingBases += codingEnd - exonStart + 1;
                        }
                    }
                    else
                    {
                        // take all of the exon's bases
                        totalCodingBases += exonEnd - exonStart + 1;

                        if (position >= exonStart)
                        {
                            if (position < exonEnd)
                                codingBases += position - exonStart + 1;
                            else
                                codingBases += exonEnd - exonStart + 1;
                        }
                    }
                }
            }

            if (!isForwardStrand)
            {
                codingBases = totalCodingBases - codingBases;
            }
        }

        Transcript transcript = new Transcript(geneAnnotation, first.TransId, first.TransName,
                upExonRank, upExonPhase, downExonRank, downExonPhase,
                (int)codingBases, (int)totalCodingBases,
                exonMax, first.IsCanonical, first.TransStart, first.TransEnd,
                first.CodingStart, first.CodingEnd);

        transcript.setBioType(first.BioType);
        transcript.setExonDistances((int)nextUpDistance, (int)nextDownDistance);

        if(isCodingTypeOverride)
            transcript.setCodingType(TRANS_CODING_TYPE_CODING);

        return transcript;
    }

    public static int EXON_RANK_MIN = 0;
    public static int EXON_RANK_MAX = 1;

    public int[] getExonRankings(final String geneId, long position)
    {
        // finds the exon before and after this position, setting to -1 if before the first or beyond the last exon
        int[] exonData = new int[EXON_RANK_MAX + 1];

        final List<TranscriptExonData> exonDataList = getTranscriptExons(geneId, "");

        if (exonDataList == null || exonDataList.isEmpty())
            return exonData;

        return getExonRankings(exonDataList, position);
    }

    public static int[] getExonRankings(final List<TranscriptExonData> exonDataList, long position)
    {
        int[] exonData = new int[EXON_RANK_MAX + 1];

        // first test a position outside the range of the exons
        final TranscriptExonData firstExon = exonDataList.get(0);
        final TranscriptExonData lastExon = exonDataList.get(exonDataList.size() - 1);

        if((position < firstExon.ExonStart && firstExon.Strand == 1) || (position > lastExon.ExonEnd && lastExon.Strand == -1))
        {
            exonData[EXON_RANK_MIN] = 0;
            exonData[EXON_RANK_MAX] = 1;
        }
        else if((position < firstExon.ExonStart && firstExon.Strand == -1) || (position > lastExon.ExonEnd && lastExon.Strand == 1))
        {
            exonData[EXON_RANK_MIN] = exonDataList.size();
            exonData[EXON_RANK_MAX] = -1;
        }
        else
        {
            for(int i = 0; i < exonDataList.size(); ++i)
            {
                final TranscriptExonData transExonData = exonDataList.get(i);
                final TranscriptExonData nextTransExonData = i < exonDataList.size() - 1 ? exonDataList.get(i+1) : null;

                if(position >= transExonData.ExonStart && position <= transExonData.ExonEnd)
                {
                    exonData[EXON_RANK_MIN] = transExonData.ExonRank;
                    exonData[EXON_RANK_MAX] = transExonData.ExonRank;
                    break;
                }

                if(nextTransExonData != null && position > transExonData.ExonEnd && position < nextTransExonData.ExonStart)
                {
                    if(transExonData.Strand == 1)
                    {
                        exonData[EXON_RANK_MIN] = transExonData.ExonRank;
                        exonData[EXON_RANK_MAX] = nextTransExonData.ExonRank;
                    }
                    else
                    {
                        exonData[EXON_RANK_MIN] = nextTransExonData.ExonRank;
                        exonData[EXON_RANK_MAX] = transExonData.ExonRank;
                    }

                    break;
                }
            }
        }

        return exonData;
    }

    public static int PSEUDO_GENE_DATA_TRANS_ID = 0;
    public static int PSEUDO_GENE_DATA_EXON_RANK = 1;
    public static int PSEUDO_GENE_DATA_EXON_MAX = 2;
    public static int PSEUDO_GENE_DATA_EXON_LENGTH = 3;

    public String[] getExonDetailsForPosition(final GeneAnnotation gene, long posStart, long posEnd)
    {
        String[] exonMatchData = new String[PSEUDO_GENE_DATA_EXON_LENGTH +1];
        exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] = null;

        List<TranscriptExonData> transExonDataList = getTransExonData(gene.StableId);

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            int exonCount = transcriptExons.size();
            for (int i = 0; i < exonCount; ++i)
            {
                final TranscriptExonData exonData = transcriptExons.get(i);

                if(abs(exonData.ExonStart - posStart) > 1 || abs(exonData.ExonEnd - posEnd) > 1)
                    continue;

                // found a match
                if (exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] == null || exonData.IsCanonical)
                {
                    exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] = exonData.TransName;
                    exonMatchData[PSEUDO_GENE_DATA_EXON_RANK] = Integer.toString(exonData.ExonRank);
                    exonMatchData[PSEUDO_GENE_DATA_EXON_MAX] = Integer.toString(exonCount);
                    exonMatchData[PSEUDO_GENE_DATA_EXON_LENGTH] = Long.toString(exonData.ExonEnd - exonData.ExonStart);

                    if (exonData.IsCanonical)
                        break;
                }
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        return exonMatchData;
    }

    public boolean loadEnsemblData(boolean delayTranscriptLoading)
    {
        if(!EnsemblDAO.loadEnsemblGeneData(mDataPath, mChromosomeGeneDataMap, mChromosomeReverseGeneDataMap))
            return false;

        if(!delayTranscriptLoading)
        {
            if (!EnsemblDAO.loadTranscriptExonData(mDataPath, mGeneTransExonDataMap, Maps.newHashMap()))
                return false;

            if(!EnsemblDAO.loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, Maps.newHashMap()))
                return false;
        }

        return true;
    }

    public boolean loadEnsemblTranscriptData(final Map<String,Boolean> uniqueGeneIds)
    {
        mGeneTransExonDataMap.clear();

        if(!EnsemblDAO.loadTranscriptExonData(mDataPath, mGeneTransExonDataMap, uniqueGeneIds))
            return false;

        Map<Integer,Boolean> uniqueTransIds = Maps.newHashMap();

        for(List<TranscriptExonData> transExonList : mGeneTransExonDataMap.values())
        {
            for(TranscriptExonData transExonData : transExonList)
            {
                uniqueTransIds.put(transExonData.TransId, true);
            }
        }

        return EnsemblDAO.loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, uniqueTransIds);
    }

    public void writeBreakendData(final String sampleId, final List<StructuralVariantAnnotation> annotations)
    {
        if(mDataPath.isEmpty())
            return;

        LOGGER.debug("writing {} breakend data to file", annotations.size());

        try
        {
            if(mBreakendWriter == null)
            {
                mBreakendWriter = createBufferedWriter(mDataPath + "SV_BREAKENDS.csv", false);

                mBreakendWriter.write("SampleId,SvId,IsStart,Chromosome,Position,Orientation,Type");
                mBreakendWriter.write(",GeneName,GeneStableId,GeneStrand,TranscriptId,IsCanonical,BioType,TransStart,TransEnd");
                mBreakendWriter.write(",ExonRankUp,ExonPhaseUp,ExonRankDown,ExonPhaseDown,CodingBases,TotalCodingBases,Disruptive");
                mBreakendWriter.write(",ExonMax,CodingStart,CodingEnd,RegionType,CodingType,ExonDistanceUp,ExonDistanceDown");
                mBreakendWriter.newLine();
            }

            BufferedWriter writer = mBreakendWriter;

            for(final StructuralVariantAnnotation annotation : annotations)
            {
                if(annotation.annotations().isEmpty())
                    continue;

                for(final GeneAnnotation geneAnnotation : annotation.annotations())
                {
                    for(final Transcript transcript : geneAnnotation.transcripts())
                    {
                        final StructuralVariant var = annotation.variant();

                        boolean isStart = geneAnnotation.isStart();

                        writer.write(String.format("%s,%d,%s,%s,%d,%d,%s",
                                sampleId, var.primaryKey(), isStart, var.chromosome(isStart), var.position(isStart),
                                var.orientation(isStart), var.type()));

                        // Gene info: geneName, geneStableId, geneStrand, transcriptId
                        writer.write(
                                String.format(",%s,%s,%d,%s,%s,%s,%d,%d",
                                        geneAnnotation.GeneName, geneAnnotation.StableId, geneAnnotation.Strand,
                                        transcript.StableId, transcript.isCanonical(), transcript.bioType(),
                                        transcript.TranscriptStart, transcript.TranscriptEnd));

                        // Transcript info: exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonStart, exonEnd, exonMax, canonical, codingStart, codingEnd
                        writer.write(
                                String.format(",%d,%d,%d,%d,%d,%d,%s",
                                        transcript.ExonUpstream, transcript.ExonUpstreamPhase,
                                        transcript.ExonDownstream, transcript.ExonDownstreamPhase,
                                        transcript.codingBases(), transcript.totalCodingBases(), transcript.isDisruptive()));

                        writer.write(
                                String.format(",%d,%d,%d,%s,%s,%d,%d",
                                        transcript.ExonMax, transcript.codingStart(), transcript.codingEnd(),
                                        transcript.regionType(), transcript.codingType(),
                                        transcript.exonDistanceUp(), transcript.exonDistanceDown()));

                        writer.newLine();
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing breakend data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mBreakendWriter);
    }

    public void writeGeneProbabilityData()
    {
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChromosomeGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData :entry.getValue())
            {
                final List<TranscriptExonData> transExonDataList = getTransExonData(geneData.GeneId);

                setGenePhasingCounts(geneData, transExonDataList);
            }
        }
    }

    public static void setGenePhasingCounts(final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList)
    {
        long geneStart = geneData.GeneStart;
        long geneEnd = geneData.GeneEnd;
        int geneLength = (int) (geneEnd - geneStart + 1);

        boolean[][] geneBases = new boolean[geneLength][GENE_PHASING_REGION_MAX];

        boolean hasCodingExons = false;
        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            if (geneData.Strand == 1)
            {
                for (int i = 0; i < transcriptExons.size(); ++i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    if (exonData.CodingStart == null)
                        break;

                    hasCodingExons = true;

                    if (exonData.ExonStart > exonData.CodingEnd) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == 0)
                    {
                        for (long j = geneStart; j <= exonData.CodingStart; ++j)
                        {
                            int gbPos = (int) (j - geneStart);
                            geneBases[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonEnd < exonData.CodingStart) // already accounted for
                        continue;

                    // now handle the exon's phasing
                    long codingStart = max(exonData.ExonStart, exonData.CodingStart);
                    for (long j = codingStart; j <= exonData.ExonEnd; ++j)
                    {
                        int gbPos = (int) (j - geneStart);

                        if (j > exonData.CodingEnd)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (j - codingStart);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i < transcriptExons.size() - 1)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i + 1);

                        if (nextExon.ExonStart <= nextExon.CodingEnd)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = exonData.ExonEnd + 1; j < nextExon.ExonStart; ++j)
                            {
                                int gbPos = (int) (j - geneStart);
                                geneBases[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }
            else
            {
                // navigate through as per the exon rank
                for (int i = transcriptExons.size() - 1; i >= 0; --i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    if (exonData.CodingStart == null)
                        break;

                    hasCodingExons = true;

                    if (exonData.ExonEnd < exonData.CodingStart) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == transcriptExons.size() - 1)
                    {
                        for (long j = exonData.CodingEnd; j <= geneEnd; ++j)
                        {
                            int gbPos = (int) (geneEnd - j);
                            geneBases[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonStart > exonData.CodingEnd) // already accounted for
                        continue;

                    // allocate the exon's phasing - working backwwards this time
                    long codingEnd = min(exonData.ExonEnd, exonData.CodingEnd);
                    for (long j = codingEnd; j >= exonData.ExonStart; --j)
                    {
                        int gbPos = (int) (geneEnd - j);

                        if (j < exonData.CodingStart)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (codingEnd - j);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i > 0)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i - 1);

                        if (nextExon.ExonEnd >= nextExon.CodingStart)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = nextExon.ExonEnd + 1; j < exonData.ExonStart; ++j)
                            {
                                int gbPos = (int) (geneEnd - j);
                                geneBases[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        // now compute the number of bases for each phasing region
        int[] regionTotals = geneData.getRegionTotals();

        for(int i = 0; i < geneLength; ++i)
        {
            for(int j = 0; j < GENE_PHASING_REGION_MAX; ++j)
            {
                if(geneBases[i][j])
                    ++regionTotals[j];
            }
        }

        if(hasCodingExons)
        {
            LOGGER.debug("gene({}) length({}) region counts: pre-coding({}) phases(0={} 1={} 2={})",
                    geneData.GeneId, geneData.GeneName, geneLength, regionTotals[GENE_PHASING_REGION_5P_UTR],
                    regionTotals[GENE_PHASING_REGION_CODING_0], regionTotals[GENE_PHASING_REGION_CODING_1],
                    regionTotals[GENE_PHASING_REGION_CODING_2]);
        }
    }


}
