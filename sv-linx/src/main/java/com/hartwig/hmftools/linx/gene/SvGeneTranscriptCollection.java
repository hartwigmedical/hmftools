package com.hartwig.hmftools.linx.gene;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.Transcript.TRANS_CODING_TYPE_CODING;
import static com.hartwig.hmftools.linx.gene.EnsemblDAO.ENSEMBL_TRANS_SPLICE_DATA_FILE;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

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
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptProteinData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<String, List<TranscriptData>> mTranscriptDataMap;
    private Map<String, List<EnsemblGeneData>> mChrGeneDataMap;
    private Map<Integer, List<TranscriptProteinData>> mEnsemblProteinDataMap;
    private Map<Integer,Long> mTransSpliceAcceptorPosDataMap;
    private Map<String, EnsemblGeneData> mGeneDataMap; // keyed by geneId
    private Map<String, EnsemblGeneData> mGeneNameIdMap; // for faster look-up by name

    // the distance upstream of a gene for a breakend to be consider a fusion candidate
    public static int PRE_GENE_PROMOTOR_DISTANCE = 100000;

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mTranscriptDataMap = Maps.newHashMap();
        mChrGeneDataMap = Maps.newHashMap();
        mEnsemblProteinDataMap = Maps.newHashMap();
        mTransSpliceAcceptorPosDataMap = Maps.newHashMap();
        mGeneDataMap = Maps.newHashMap();
        mGeneNameIdMap = Maps.newHashMap();
    }

    public void setDataPath(final String dataPath)
    {
        mDataPath = dataPath;

        if(!mDataPath.endsWith(File.separator))
            mDataPath += File.separator;
    }

    public final Map<String, List<TranscriptData>> getTranscriptDataMap() { return mTranscriptDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrGeneDataMap() { return mChrGeneDataMap; }
    public Map<Integer, List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }
    public Map<Integer,Long> getTransSpliceAcceptorPosDataMap() { return mTransSpliceAcceptorPosDataMap; }

    public final EnsemblGeneData getGeneDataByName(final String geneName)
    {
        if(!mGeneNameIdMap.isEmpty())
            return mGeneNameIdMap.get(geneName);

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
        if(!mGeneDataMap.isEmpty())
            return;

        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                mGeneDataMap.put(geneData.GeneId, geneData);
            }
        }
    }

    public void createGeneNameIdMap()
    {
        if(!mGeneNameIdMap.isEmpty())
            return;

        for(Map.Entry<String, List<EnsemblGeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                mGeneNameIdMap.put(geneData.GeneName, geneData);
            }
        }
    }

    public List<TranscriptData> getTranscripts(final String geneId)
    {
        return mTranscriptDataMap.get(geneId);
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

        if(geneRegions == null)
            return geneAnnotations;

        final List<EnsemblGeneData> matchedGenes = findGeneRegions(position, geneRegions, upstreamDistance);

        // now look up relevant transcript and exon information
        for(final EnsemblGeneData geneData : matchedGenes)
        {
            final List<TranscriptData> transcriptDataList = mTranscriptDataMap.get(geneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            GeneAnnotation currentGene = new GeneAnnotation(svId, isStart, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.Synonyms, geneData.KaryotypeBand);

            currentGene.setGeneData(geneData);

            // collect up all the relevant exons for each unique transcript to analyse as a collection

            for(TranscriptData transData : transcriptDataList)
            {
                Transcript transcript = extractTranscriptExonData(transData, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);

                    setAlternativeTranscriptPhasings(transcript, transData.exons(), position, orientation);

                    // annotate with preceding gene info if the up distance isn't set
                    if(transcript.exonDistanceUp() == -1)
                    {
                        setPrecedingGeneDistance(transcript, position);
                    }
                }
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    private void setPrecedingGeneDistance(Transcript transcript, long position)
    {
        // annotate with preceding gene info if the up distance isn't set
        long precedingGeneSAPos = findPrecedingGeneSpliceAcceptorPosition(transcript.TransId);

        if(precedingGeneSAPos >= 0)
        {
            long preDistance = transcript.parent().Strand == 1 ? position - precedingGeneSAPos : precedingGeneSAPos - position;
            transcript.setExonDistances((int)preDistance, transcript.exonDistanceDown());
        }
    }

    public final TranscriptData getTranscriptData(final String geneId, final String transcriptId)
    {
        final List<TranscriptData> transDataList = mTranscriptDataMap.get(geneId);

        if (transDataList == null || transDataList.isEmpty())
            return null;

        for(final TranscriptData transData : transDataList)
        {
            if(transcriptId.isEmpty() && transData.IsCanonical)
                return transData;
            else if(transData.TransName.equals(transcriptId))
                return transData;
        }

        return null;
    }

    public final List<EnsemblGeneData> findGenesByRegion(final String chromosome, long posStart, long posEnd)
    {
        // find genes if any of their transcripts are within this position
        List<EnsemblGeneData> genesList = Lists.newArrayList();

        final List<EnsemblGeneData> geneDataList = mChrGeneDataMap.get(chromosome);

        for(final EnsemblGeneData geneData : geneDataList)
        {
            if(posStart > geneData.GeneEnd || posEnd < geneData.GeneStart)
                continue;

            final List<TranscriptData> transList = mTranscriptDataMap.get(geneData.GeneId);

            if(transList == null || transList.isEmpty())
                continue;

            for(final TranscriptData transData : transList)
            {
                if (posStart <= transData.TransStart && posEnd >= transData.TransEnd)
                {
                    genesList.add(geneData);
                    break;
                }
            }
        }

        return genesList;
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

                TranscriptData transData = getTranscriptData(gene.GeneId, "");

                if(transData == null || transData.exons().size() <= 1)
                    continue;

                // otherwise taken the start of the last exon
                ExonData lastExonData = transData.exons().get(transData.exons().size() - 1);
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

                TranscriptData transData = getTranscriptData(gene.GeneId, "");

                if(transData == null || transData.exons().size() <= 1)
                    continue;

                // otherwise taken the start of the last exon
                ExonData lastExonData = transData.exons().get(0);
                return lastExonData.ExonEnd;
            }
        }

        return -1;
    }

    public static Transcript extractTranscriptExonData(final TranscriptData transData, long position,
            final GeneAnnotation geneAnnotation)
    {
        final List<ExonData> exonList = transData.exons();
        int exonMax = exonList.size();

        boolean isForwardStrand = geneAnnotation.Strand == 1;

        int upExonRank = -1;
        int upExonPhase = -1;
        int downExonRank = -1;
        int downExonPhase = -1;
        long nextUpDistance = -1;
        long nextDownDistance = -1;
        boolean isCodingTypeOverride = false;

        // first check for a position outside the exon boundaries
        final ExonData firstExon = exonList.get(0);
        final ExonData lastExon = exonList.get(exonList.size()-1);

        // for forward-strand transcripts the current exon is downstream, the previous is upstream
        // and the end-phase is taken from the upstream previous exon, the phase from the current downstream exon

        // for reverse-strand transcripts the current exon is upstream, the previous is downstream
        // and the end-phase is taken from the upstream (current) exon, the phase from the downstream (previous) exon

        // for each exon, the 'phase' is always the phase at the start of the exon in the direction of transcrition
        // regardless of strand direction, and 'end_phase' is the phase at the end of the exon

        if(position < firstExon.ExonStart)
        {
            if(isForwardStrand)
            {
                // proceed to the next exon assuming its splice acceptor is required
                final ExonData firstSpaExon = exonList.size() > 1 ? exonList.get(1) : firstExon;
                downExonRank = firstSpaExon.ExonRank;
                downExonPhase = firstSpaExon.ExonPhase;
                nextDownDistance = firstSpaExon.ExonStart - position;

                // correct the phasing if the next exon starts the coding region
                if(transData.CodingStart != null && firstSpaExon.ExonStart == transData.CodingStart)
                    downExonPhase = -1;

                isCodingTypeOverride = transData.CodingStart != null && firstSpaExon.ExonStart > transData.CodingStart;

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
                final ExonData firstSpaExon = exonList.size() > 1 ? exonList.get(exonList.size()-2) : lastExon;
                downExonRank = firstSpaExon.ExonRank;
                downExonPhase = firstSpaExon.ExonPhase;
                nextDownDistance = position - lastExon.ExonEnd;

                if(transData.CodingEnd != null && firstSpaExon.ExonEnd == transData.CodingEnd)
                    downExonPhase = -1;

                isCodingTypeOverride = transData.CodingEnd != null && firstSpaExon.ExonEnd < transData.CodingEnd;

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
            for (int index = 0; index < exonList.size(); ++index)
            {
                final ExonData exonData = exonList.get(index);

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

                        if(index < exonList.size() - 1)
                        {
                            final ExonData nextExonData = exonList.get(index + 1);
                            nextDownDistance = nextExonData.ExonStart - position;
                        }
                    }
                    else
                    {
                        nextUpDistance = exonData.ExonEnd - position;

                        if(index > 1)
                        {
                            // first splice acceptor is the second exon (or later on)
                            final ExonData prevExonData = exonList.get(index - 1);
                            nextDownDistance = position - prevExonData.ExonEnd;
                        }
                    }

                    break;
                }
                else if(position < exonData.ExonStart)
                {
                    // position falls between this exon and the previous one
                    final ExonData prevExonData = exonList.get(index-1);

                    if(isForwardStrand)
                    {
                        // the current exon is downstream, the prevous one is upstream
                        upExonRank = prevExonData.ExonRank;
                        upExonPhase = prevExonData.ExonPhaseEnd;
                        downExonRank = exonData.ExonRank;
                        downExonPhase = exonData.ExonPhase;
                        nextDownDistance = exonData.ExonStart - position;
                        nextUpDistance = position - prevExonData.ExonEnd;

                        if(transData.CodingStart != null && exonData.ExonStart == transData.CodingStart)
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

                        if(transData.CodingEnd != null && prevExonData.ExonEnd == transData.CodingEnd)
                            downExonPhase = -1;
                    }

                    break;
                }
            }
        }

        // now calculate coding bases for this transcript
        // for the given position, determine how many coding bases occur prior to the position
        // in the direction of the transcript

        boolean isCoding = transData.CodingStart != null && transData.CodingEnd != null;
        long codingStart = transData.CodingStart != null ? transData.CodingStart : 0;
        long codingEnd = transData.CodingEnd != null ? transData.CodingEnd : 0;
        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        long codingBases = 0;
        long totalCodingBases = 0;

        if(isCoding)
        {
            for (ExonData exonData : exonList)
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

        Transcript transcript = new Transcript(geneAnnotation, transData.TransId, transData.TransName,
                upExonRank, upExonPhase, downExonRank, downExonPhase,
                (int)codingBases, (int)totalCodingBases,
                exonMax, transData.IsCanonical, transData.TransStart, transData.TransEnd,
                transData.CodingStart, transData.CodingEnd);

        transcript.setBioType(transData.BioType);
        transcript.setExonDistances((int)nextUpDistance, (int)nextDownDistance);

        if(isCodingTypeOverride)
            transcript.setCodingType(TRANS_CODING_TYPE_CODING);

        return transcript;
    }

    public static int EXON_RANK_MIN = 0;
    public static int EXON_RANK_MAX = 1;
    public static int EXON_PHASE_MIN = 2;
    public static int EXON_PHASE_MAX = 3;

    public int[] getExonRankings(final String geneId, long position)
    {
        // finds the exon before and after this position, setting to -1 if before the first or beyond the last exon
        int[] exonData = new int[EXON_PHASE_MAX + 1];

        final TranscriptData transData = getTranscriptData(geneId, "");

        if (transData == null || transData.exons().isEmpty())
            return exonData;

        return getExonRankings(transData.Strand, transData.exons(), position);
    }

    public static int[] getExonRankings(int strand, final List<ExonData> exonDataList, long position)
    {
        int[] exonData = new int[EXON_PHASE_MAX + 1];

        // first test a position outside the range of the exons
        final ExonData firstExon = exonDataList.get(0);
        final ExonData lastExon = exonDataList.get(exonDataList.size() - 1);

        if((position < firstExon.ExonStart && strand == 1) || (position > lastExon.ExonEnd && strand == -1))
        {
            // before the start of the transcript
            exonData[EXON_RANK_MIN] = 0;
            exonData[EXON_RANK_MAX] = 1;
            exonData[EXON_PHASE_MIN] = -1;
            exonData[EXON_PHASE_MAX] = -1;
        }
        else if((position < firstExon.ExonStart && strand == -1) || (position > lastExon.ExonEnd && strand == 1))
        {
            // past the end of the transcript
            exonData[EXON_RANK_MIN] = exonDataList.size();
            exonData[EXON_RANK_MAX] = -1;
            exonData[EXON_PHASE_MIN] = -1;
            exonData[EXON_PHASE_MAX] = -1;
        }
        else
        {
            for(int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData transExonData = exonDataList.get(i);
                final ExonData nextTransExonData = i < exonDataList.size() - 1 ? exonDataList.get(i+1) : null;

                if(position == transExonData.ExonEnd || position == transExonData.ExonStart)
                {
                    // position matches the bounds of an exon
                    exonData[EXON_RANK_MIN] = transExonData.ExonRank;
                    exonData[EXON_RANK_MAX] = transExonData.ExonRank;

                    if((strand == 1) == (position == transExonData.ExonStart))
                    {
                        exonData[EXON_PHASE_MIN] = transExonData.ExonPhase;
                        exonData[EXON_PHASE_MAX] = transExonData.ExonPhase;
                    }
                    else
                    {
                        exonData[EXON_PHASE_MIN] = transExonData.ExonPhaseEnd;
                        exonData[EXON_PHASE_MAX] = transExonData.ExonPhaseEnd;
                    }
                    break;
                }

                if(position >= transExonData.ExonStart && position <= transExonData.ExonEnd)
                {
                    // position matches within or at the bounds of an exon
                    exonData[EXON_RANK_MIN] = transExonData.ExonRank;
                    exonData[EXON_RANK_MAX] = transExonData.ExonRank;
                    exonData[EXON_PHASE_MIN] = transExonData.ExonPhase;
                    exonData[EXON_PHASE_MAX] = transExonData.ExonPhase;
                    break;
                }

                if(nextTransExonData != null && position > transExonData.ExonEnd && position < nextTransExonData.ExonStart)
                {
                    if(strand == 1)
                    {
                        exonData[EXON_RANK_MIN] = transExonData.ExonRank;
                        exonData[EXON_RANK_MAX] = nextTransExonData.ExonRank;
                        exonData[EXON_PHASE_MIN] = transExonData.ExonPhase;
                        exonData[EXON_PHASE_MAX] = nextTransExonData.ExonPhase;
                    }
                    else
                    {
                        exonData[EXON_RANK_MIN] = nextTransExonData.ExonRank;
                        exonData[EXON_RANK_MAX] = transExonData.ExonRank;
                        exonData[EXON_PHASE_MIN] = nextTransExonData.ExonPhase;
                        exonData[EXON_PHASE_MAX] = transExonData.ExonPhase;
                    }

                    break;
                }
            }
        }

        return exonData;
    }

    public static void setAlternativeTranscriptPhasings(Transcript transcript, final List<ExonData> exonDataList,
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
            final ExonData exonData = exonDataList.get(i);

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

    public static int PSEUDO_GENE_DATA_TRANS_ID = 0;
    public static int PSEUDO_GENE_DATA_EXON_RANK = 1;
    public static int PSEUDO_GENE_DATA_EXON_MAX = 2;
    public static int PSEUDO_GENE_DATA_EXON_LENGTH = 3;

    public String[] getExonDetailsForPosition(final GeneAnnotation gene, long posStart, long posEnd)
    {
        String[] exonMatchData = new String[PSEUDO_GENE_DATA_EXON_LENGTH +1];
        exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] = null;

        List<TranscriptData> transDataList = getTranscripts(gene.StableId);

        for(final TranscriptData transData : transDataList)
        {
            int exonCount = transData.exons().size();

            for (int i = 0; i < exonCount; ++i)
            {
                final ExonData exonData = transData.exons().get(i);

                if(abs(exonData.ExonStart - posStart) > 1 || abs(exonData.ExonEnd - posEnd) > 1)
                    continue;

                // found a match
                if (exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] == null || transData.IsCanonical)
                {
                    exonMatchData[PSEUDO_GENE_DATA_TRANS_ID] = transData.TransName;
                    exonMatchData[PSEUDO_GENE_DATA_EXON_RANK] = Integer.toString(exonData.ExonRank);
                    exonMatchData[PSEUDO_GENE_DATA_EXON_MAX] = Integer.toString(exonCount);
                    exonMatchData[PSEUDO_GENE_DATA_EXON_LENGTH] = Long.toString(exonData.ExonEnd - exonData.ExonStart);

                    if (transData.IsCanonical)
                        break;
                }
            }
        }

        return exonMatchData;
    }

    public boolean loadEnsemblData(boolean delayTranscriptLoading)
    {
        if(!EnsemblDAO.loadEnsemblGeneData(mDataPath, mChrGeneDataMap))
            return false;

        if(!delayTranscriptLoading)
        {
            if(!EnsemblDAO.loadTranscriptData(mDataPath, mTranscriptDataMap, Lists.newArrayList()))
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
                // createTranscriptPreGenePositionData(getChrGeneDataMap(), getGeneExonDataMap());
            }
        }

        return true;
    }

    public boolean loadEnsemblTranscriptData(final List<String> uniqueGeneIds)
    {
        if(!EnsemblDAO.loadTranscriptData(mDataPath, mTranscriptDataMap, uniqueGeneIds))
            return false;

        List<Integer> uniqueTransIds = Lists.newArrayList();

        for(List<TranscriptData> transDataList : mTranscriptDataMap.values())
        {
            for(TranscriptData transData : transDataList)
            {
                if(!uniqueTransIds.contains(transData.TransId))
                    uniqueTransIds.add(transData.TransId);
            }
        }

        if(!EnsemblDAO.loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, uniqueTransIds))
            return false;

        return EnsemblDAO.loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, uniqueTransIds);
    }

    public void createTranscriptPreGenePositionData(final Map<String, List<EnsemblGeneData>> chrGeneDataMap,
            final Map<String, List<TranscriptData>> transcriptDataMap)
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
                    final List<TranscriptData> transDataList = transcriptDataMap.get(gene.GeneId);

                    if (transDataList == null || transDataList.isEmpty())
                        continue;


                    for(final TranscriptData transData : transDataList)
                    {
                        long transStartPos = gene.Strand == 1 ? transData.TransStart : transData.TransEnd;

                        long firstSpliceAcceptorPos = findFirstSpliceAcceptor(transStartPos, gene.Strand, proximateGenes, transcriptDataMap);

                        long distance = -1;

                        if(firstSpliceAcceptorPos >= 0)
                            distance = gene.Strand == 1 ? transStartPos - firstSpliceAcceptorPos : firstSpliceAcceptorPos - transStartPos;

                        // cache value and continue
                        writer.write(String.format("%s,%d,%s,%d,%d,%d",
                                gene.GeneId, transData.TransId, transData.TransName, transStartPos, firstSpliceAcceptorPos, distance));

                        writer.newLine();

                        mTransSpliceAcceptorPosDataMap.put(transData.TransId, firstSpliceAcceptorPos);
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
            long transStartPos, int strand, final List<String> proximateGenes, final Map<String, List<TranscriptData>> transDataMap)
    {
        long closestPosition = -1;

        for(final String geneId : proximateGenes)
        {
            final List<TranscriptData> transDataList = transDataMap.get(geneId);

            if(transDataList.isEmpty())
                continue;

            for(final TranscriptData transData : transDataList)
            {
                for (final ExonData exonData : transData.exons())
                {
                    // check if exon is upstream and if so record its position
                    if (strand == 1 && exonData.ExonEnd < transStartPos)
                    {
                        if (closestPosition == -1 || exonData.ExonStart > closestPosition)
                            closestPosition = exonData.ExonStart;
                    }
                    else if (strand == -1 && exonData.ExonStart > transStartPos)
                    {
                        if (closestPosition == -1 || exonData.ExonEnd < closestPosition)
                            closestPosition = exonData.ExonEnd;
                    }
                }
            }
        }

        return closestPosition;
    }

    public static Long[] getProteinDomainPositions(final TranscriptProteinData proteinData, final TranscriptData transData)
    {
        Long[] domainPositions = {null, null};

        if(transData.exons().isEmpty())
            return domainPositions;

        Long codingStart = transData.CodingStart;
        Long codingEnd = transData.CodingEnd;

        if(codingStart == null || codingEnd == null)
            return domainPositions;

        int preProteinBases = proteinData.SeqStart * 3;
        int proteinBases = (proteinData.SeqEnd - proteinData.SeqStart) * 3;

        long proteinStart = -1;
        long proteinEnd = -1;

        if(transData.Strand == 1)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exonData = transData.exons().get(i);

                if(exonData.ExonEnd < codingStart)
                    continue;

                if(preProteinBases > 0)
                {
                    long refStartPos = max(codingStart, exonData.ExonStart);
                    long exonCodingBases = exonData.ExonEnd - refStartPos;

                    if(exonCodingBases >= preProteinBases)
                    {
                        proteinStart = refStartPos + preProteinBases;
                        preProteinBases = 0;
                    }
                    else
                    {
                        preProteinBases -= exonCodingBases;
                        continue;
                    }
                }

                long startPos = max(exonData.ExonStart, proteinStart);
                long exonBases = exonData.ExonEnd - startPos;

                if(exonBases >= proteinBases)
                {
                    proteinEnd = startPos + proteinBases;
                    break;
                }
                else
                {
                    proteinBases -= exonBases;
                }
            }
        }
        else
        {
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                final ExonData exonData = transData.exons().get(i);

                if(exonData.ExonStart > codingEnd)
                    continue;

                if(preProteinBases > 0)
                {
                    long refStartPos = min(codingEnd, exonData.ExonEnd);
                    long exonCodingBases = refStartPos - exonData.ExonStart;

                    if(exonCodingBases >= preProteinBases)
                    {
                        proteinEnd = refStartPos - preProteinBases;
                        preProteinBases = 0;
                    }
                    else
                    {
                        preProteinBases -= exonCodingBases;
                        continue;
                    }
                }

                long startPos = min(exonData.ExonEnd, proteinEnd);
                long exonBases = startPos - exonData.ExonStart;

                if(exonBases >= proteinBases)
                {
                    proteinStart = startPos - proteinBases;
                    break;
                }
                else
                {
                    proteinBases -= exonBases;
                }
            }
        }

        if(proteinEnd == -1 || proteinStart == -1)
            return domainPositions;

        domainPositions[SE_START] = proteinStart;
        domainPositions[SE_END] = proteinEnd;

        return domainPositions;
    }

}
