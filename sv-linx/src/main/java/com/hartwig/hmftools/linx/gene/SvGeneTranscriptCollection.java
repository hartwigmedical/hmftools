package com.hartwig.hmftools.linx.gene;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.TRANS_CODING_TYPE_CODING;
import static com.hartwig.hmftools.linx.gene.EnsemblDAO.ENSEMBL_TRANS_SPLICE_DATA_FILE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<String, List<TranscriptExonData>> mGeneTransExonDataMap;// keyed by GeneId (aka StableId)
    private Map<String, List<EnsemblGeneData>> mChrGeneDataMap;
    private Map<String, List<EnsemblGeneData>> mChrReverseGeneDataMap; // order by gene end instead of start, for traversal in the reverse direction
    private Map<Integer, List<TranscriptProteinData>> mEnsemblProteinDataMap;
    private Map<Integer,Long> mTransSpliceAcceptorPosDataMap;
    private Map<String, EnsemblGeneData> mGeneDataMap; // keyed by geneId

    // the distance upstream of a gene for a breakend to be consider a fusion candidate
    public static int PRE_GENE_PROMOTOR_DISTANCE = 100000;

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mGeneTransExonDataMap = Maps.newHashMap();
        mChrGeneDataMap = Maps.newHashMap();
        mChrReverseGeneDataMap = Maps.newHashMap();
        mEnsemblProteinDataMap = Maps.newHashMap();
        mTransSpliceAcceptorPosDataMap = Maps.newHashMap();
        mGeneDataMap = Maps.newHashMap();
    }

    public void setDataPath(final String dataPath)
    {
        mDataPath = dataPath;

        if(!mDataPath.endsWith(File.separator))
            mDataPath += File.separator;
    }

    public boolean hasCachedEnsemblData()
    {
        return !mChrGeneDataMap.isEmpty() && !mGeneTransExonDataMap.isEmpty();
    }

    public final Map<String, List<TranscriptExonData>> getGeneExonDataMap() { return mGeneTransExonDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrGeneDataMap() { return mChrGeneDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrReverseGeneDataMap() { return mChrReverseGeneDataMap; }
    public Map<Integer, List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }
    public Map<Integer,Long> getTransSpliceAcceptorPosDataMap() { return mTransSpliceAcceptorPosDataMap; }

    public final EnsemblGeneData getGeneDataByName(final String geneName)
    {
        return getGeneData(geneName, true);
    }

    public final EnsemblGeneData getGeneDataById(final String geneId)
    {
        if(!mGeneDataMap.isEmpty())
            return mGeneDataMap.get(geneId);

        return getGeneData(geneId, false);
    }

    private final EnsemblGeneData getGeneData(final String gene, boolean byName)
    {
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                if((byName && geneData.GeneName.equals(gene)) || (!byName && geneData.GeneId.equals(gene)))
                    return geneData;
            }
        }

        return null;
    }

    public void createGeneIdDataMap()
    {
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                mGeneDataMap.put(geneData.GeneId, geneData);
            }
        }
    }

    public List<TranscriptExonData> getTransExonData(final String geneId)
    {
        return mGeneTransExonDataMap.get(geneId);
    }

    public void populateGeneIdList(List<String> uniqueGeneIds, final String chromosome, long position, int upstreamDistance)
    {
        // find the unique set of geneIds
        final List<EnsemblGeneData> geneRegions = mChrGeneDataMap.get(chromosome);

        if (geneRegions == null)
            return;

        List<EnsemblGeneData> matchedGenes = findGeneRegions(position, geneRegions, upstreamDistance);

        for (final EnsemblGeneData geneData : matchedGenes)
        {
            if(!uniqueGeneIds.contains(geneData.GeneId))
                uniqueGeneIds.add(geneData.GeneId);
        }
    }

    public List<GeneAnnotation> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, long position,
            byte orientation, int upstreamDistance)
    {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        final List<EnsemblGeneData> geneRegions = mChrGeneDataMap.get(chromosome);
        final List<EnsemblGeneData> geneRegionsReversed = mChrReverseGeneDataMap.get(chromosome);

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

                    setAlternativeTranscriptPhasings(transcript, transcriptExons, position, orientation);

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

        final List<EnsemblGeneData> geneDataList = mChrGeneDataMap.get(chromosome);

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

    public long findPrecedingGeneSpliceAcceptorPosition(int transId)
    {
        if(mTransSpliceAcceptorPosDataMap.isEmpty())
            return -1;

        Long spliceAcceptorPos = mTransSpliceAcceptorPosDataMap.get(transId);
        return spliceAcceptorPos != null ? spliceAcceptorPos : -1;
    }

    public long findPrecedingGeneSpliceAcceptorPosition(final EnsemblGeneData refGene, List<EnsemblGeneData> geneDataList)
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

    public static Transcript extractTranscriptExonData(final List<TranscriptExonData> transcriptExons, long position,
            final GeneAnnotation geneAnnotation)
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

        // for each exon, the 'phase' is always the phase at the start of the exon regardless of strand direction,
        // and 'end_phase' is the phase at the end of the exon, but note that 'phase' will correspond to an exon end

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

                    // set distance to next and previous splice acceptor
                    if(isForwardStrand)
                    {
                        nextUpDistance = position - exonData.ExonStart;

                        if(index < transcriptExons.size() - 1)
                        {
                            final TranscriptExonData nextExonData = transcriptExons.get(index + 1);
                            nextDownDistance = nextExonData.ExonStart - position;
                        }
                    }
                    else
                    {
                        nextUpDistance = exonData.ExonEnd - position;

                        if(index > 1)
                        {
                            // first splice acceptor is the second exon (or later on)
                            final TranscriptExonData prevExonData = transcriptExons.get(index - 1);
                            nextDownDistance = position - prevExonData.ExonEnd;
                        }
                    }

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
                        // the current exon is earlier in rank
                        // the previous exon in the list has the higher rank and is dowstream
                        // the start of the next exon (ie previous here) uses 'phase' for the downstream as normal
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

    public static void setAlternativeTranscriptPhasings(Transcript transcript, final List<TranscriptExonData> exonDataList,
            long position, byte orientation)
    {
        // collect exon phasings before the position on the upstream and after it on the downstream
        boolean isUpstream = (transcript.parent().Strand * orientation) > 0;
        boolean forwardStrand = (transcript.parent().Strand == 1);

        Map<Integer,Integer> alternativePhasing = Maps.newHashMap();

        int transPhase = isUpstream ? transcript.ExonUpstreamPhase : transcript.ExonDownstreamPhase;
        int transRank = isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream;

        for(int i = 0; i < exonDataList.size(); ++i)
        {
            final TranscriptExonData exonData = exonDataList.get(i);

            if(isUpstream == forwardStrand)
            {
                if (exonData.ExonStart > position || transRank == exonData.ExonRank)
                {
                    break;
                }
            }
            else
            {
                if (position > exonData.ExonEnd || transRank == exonData.ExonRank)
                {
                    continue;
                }
            }

            int exonPhase = isUpstream ? exonData.ExonPhaseEnd : exonData.ExonPhase;
            int exonsSkipped = 0;

            if(isUpstream)
            {
                exonsSkipped = max(transRank - exonData.ExonRank, 0);
            }
            else
            {
                exonsSkipped = max(exonData.ExonRank - transRank, 0);
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

        transcript.setAlternativePhasing(alternativePhasing);
    }

    public List<Integer> getTranscriptIdsMatchingPosition(final String geneId, long position)
    {
        List<Integer> transIdList = Lists.newArrayList();

        List<TranscriptExonData> transExonDataList = getTransExonData(geneId);

        if(transExonDataList == null || transExonDataList.isEmpty())
            return transIdList;

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while(!transcriptExons.isEmpty())
        {
            for(TranscriptExonData exonData : transcriptExons)
            {
                if(abs(exonData.ExonStart - position) <= 1 || abs(exonData.ExonEnd - position) <= 1)
                {
                    transIdList.add(exonData.TransId);
                    break;
                }
            }

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        return transIdList;
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
        if(!EnsemblDAO.loadEnsemblGeneData(mDataPath, mChrGeneDataMap, mChrReverseGeneDataMap))
            return false;

        if(!delayTranscriptLoading)
        {
            if (!EnsemblDAO.loadTranscriptExonData(mDataPath, mGeneTransExonDataMap, Lists.newArrayList()))
                return false;

            if(!EnsemblDAO.loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, Lists.newArrayList()))
                return false;

            final String transSpliceFile = mDataPath + ENSEMBL_TRANS_SPLICE_DATA_FILE;

            if (Files.exists(Paths.get(transSpliceFile)))
            {
                if(!EnsemblDAO.loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, Lists.newArrayList()))
                    return false;
            }
            else
            {
                createTranscriptPreGenePositionData(getChrGeneDataMap(), getGeneExonDataMap());
            }
        }

        return true;
    }

    public boolean loadEnsemblTranscriptData(final List<String> uniqueGeneIds)
    {
        mGeneTransExonDataMap.clear();

        if(!EnsemblDAO.loadTranscriptExonData(mDataPath, mGeneTransExonDataMap, uniqueGeneIds))
            return false;

        List<Integer> uniqueTransIds = Lists.newArrayList();

        for(List<TranscriptExonData> transExonList : mGeneTransExonDataMap.values())
        {
            for(TranscriptExonData transExonData : transExonList)
            {
                if(!uniqueTransIds.contains(transExonData.TransId))
                    uniqueTransIds.add(transExonData.TransId);
            }
        }

        if(!EnsemblDAO.loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, uniqueTransIds))
            return false;

        return EnsemblDAO.loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, uniqueTransIds);
    }

    public void createTranscriptPreGenePositionData(final Map<String, List<EnsemblGeneData>> chrGeneDataMap, final Map<String, List<TranscriptExonData>> geneTransExonDataMap)
    {
        try
        {
            final String outputFile = mDataPath + ENSEMBL_TRANS_SPLICE_DATA_FILE;

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,TransId,TransName,TransStartPos,PreSpliceAcceptorPosition,Distance");
            writer.newLine();

            // for each gene, collect up any gene which overlaps it or is within the specified distance upstream from it
            for (Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
            {
                final String chromosome = entry.getKey();

                LOGGER.info("calculating pre-gene positions for chromosome({})", chromosome);

                final List<EnsemblGeneData> geneList = entry.getValue();

                for (final EnsemblGeneData gene : geneList)
                {
                    List<String> proximateGenes = Lists.newArrayList();

                    for (final EnsemblGeneData otherGene : geneList)
                    {
                        if (otherGene.Strand != gene.Strand || otherGene.GeneId.equals(gene.GeneId))
                            continue;

                        if (gene.Strand == 1)
                        {
                            if (otherGene.GeneStart >= gene.GeneStart || otherGene.GeneEnd < gene.GeneStart - PRE_GENE_PROMOTOR_DISTANCE)
                                continue;

                            proximateGenes.add(otherGene.GeneId);
                        }
                        else
                        {
                            if (otherGene.GeneEnd <= gene.GeneEnd || otherGene.GeneStart > gene.GeneEnd + PRE_GENE_PROMOTOR_DISTANCE)
                                continue;

                            proximateGenes.add(otherGene.GeneId);
                        }
                    }

                    if(proximateGenes.isEmpty())
                        continue;

                    // now set the preceding splice acceptor position for each transcript in this gene
                    final List<TranscriptExonData> transExonDataList = geneTransExonDataMap.get(gene.GeneId);

                    if (transExonDataList == null || transExonDataList.isEmpty())
                        continue;

                    int teIndex = 0;
                    List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

                    while (!transcriptExons.isEmpty())
                    {
                        final TranscriptExonData exonData = transcriptExons.get(0);
                        long transStartPos = gene.Strand == 1 ? exonData.TransStart : exonData.TransEnd;

                        long firstSpliceAcceptorPos =
                                findFirstSpliceAcceptor(transStartPos, gene.Strand, proximateGenes, geneTransExonDataMap);

                        long distance = gene.Strand == 1 ? transStartPos - firstSpliceAcceptorPos : firstSpliceAcceptorPos - transStartPos;

                        // cache value and continue
                        writer.write(String.format("%s,%d,%s,%d,%d,%d",
                                gene.GeneId, exonData.TransId, exonData.TransName, transStartPos, firstSpliceAcceptorPos, distance));

                        writer.newLine();

                        mTransSpliceAcceptorPosDataMap.put(exonData.TransId, firstSpliceAcceptorPos);

                        teIndex += transcriptExons.size();
                        transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
                    }
                }
            }

            LOGGER.info("pre-gene positions written to file: {}", outputFile);
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LOGGER.error("error writing Ensembl trans splice data file: {}", e.toString());
        }
    }

    private long findFirstSpliceAcceptor(
            long transStartPos, int strand, final List<String> proximateGenes, final Map<String, List<TranscriptExonData>> geneTransExonDataMap)
    {
        long closestPosition = -1;

        for(final String geneId : proximateGenes)
        {
            final List<TranscriptExonData> transExonDataList = geneTransExonDataMap.get(geneId);

            if(transExonDataList.isEmpty())
                continue;

            for(final TranscriptExonData exonData : transExonDataList)
            {
                // check if exon is upstream and if so record its position
                if(strand == 1 && exonData.ExonEnd < transStartPos)
                {
                    if(closestPosition == -1 || exonData.ExonStart > closestPosition)
                        closestPosition = exonData.ExonStart;
                }
                else if(strand == -1 && exonData.ExonStart > transStartPos)
                {
                    if(closestPosition == -1 || exonData.ExonEnd < closestPosition)
                        closestPosition = exonData.ExonEnd;
                }
            }
        }

        return closestPosition;
    }

}
