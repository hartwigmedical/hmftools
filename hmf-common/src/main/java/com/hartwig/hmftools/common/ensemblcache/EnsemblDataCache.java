package com.hartwig.hmftools.common.ensemblcache;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptProteinData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptSpliceAcceptorData;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.Transcript.CODING_BASES;
import static com.hartwig.hmftools.common.fusion.Transcript.TOTAL_CODING_BASES;
import static com.hartwig.hmftools.common.fusion.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class EnsemblDataCache
{
    private final String mDataPath;
    private final RefGenomeVersion mRefGenomeVersion;

    private final Map<String, List<TranscriptData>> mTranscriptDataMap;
    private final Map<String, List<EnsemblGeneData>> mChrGeneDataMap;
    private final Map<Integer, List<TranscriptProteinData>> mEnsemblProteinDataMap;
    private final Map<Integer,Integer> mTransSpliceAcceptorPosDataMap;
    private final Map<String,EnsemblGeneData> mGeneDataMap; // keyed by geneId
    private final Map<String,EnsemblGeneData> mGeneNameIdMap; // for faster look-up by name

    // whether to load more details information for each transcript - exons, protein domains, splice positions etc
    private boolean mRequireExons;
    private boolean mRequireProteinDomains;
    private boolean mRequireSplicePositions;
    private boolean mCanonicalTranscriptsOnly;

    private final Map<EnsemblGeneData,Integer> mDownstreamGeneAnnotations;

    private final List<String> mRestrictedGeneIdList = Lists.newArrayList();

    public EnsemblDataCache(final String dataPath, final RefGenomeVersion refGenomeVersion)
    {
        mDataPath = dataPath.endsWith(File.separator) ? dataPath : dataPath + File.separator;
        mRefGenomeVersion = refGenomeVersion;

        mTranscriptDataMap = Maps.newHashMap();
        mChrGeneDataMap = Maps.newHashMap();
        mEnsemblProteinDataMap = Maps.newHashMap();
        mTransSpliceAcceptorPosDataMap = Maps.newHashMap();
        mGeneDataMap = Maps.newHashMap();
        mGeneNameIdMap = Maps.newHashMap();
        mRequireExons = true;
        mRequireProteinDomains = false;
        mRequireSplicePositions = false;
        mCanonicalTranscriptsOnly = false;
        mDownstreamGeneAnnotations = Maps.newHashMap();
    }

    public void setRestrictedGeneIdList(final List<String> geneIds)
    {
        mRestrictedGeneIdList.clear();
        mRestrictedGeneIdList.addAll(geneIds);
    }

    public final Map<EnsemblGeneData,Integer> getDownstreamGeneAnnotations() { return mDownstreamGeneAnnotations; }

    public void setRequiredData(boolean exons, boolean proteinDomains, boolean splicePositions, boolean canonicalOnly)
    {
        mRequireExons = exons;
        mRequireSplicePositions = splicePositions;
        mRequireProteinDomains = proteinDomains;
        mCanonicalTranscriptsOnly = canonicalOnly;
    }

    public final Map<String, List<TranscriptData>> getTranscriptDataMap() { return mTranscriptDataMap; }
    public final Map<String, List<EnsemblGeneData>> getChrGeneDataMap() { return mChrGeneDataMap; }
    public Map<Integer, List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }
    public Map<Integer,Integer> getTransSpliceAcceptorPosDataMap() { return mTransSpliceAcceptorPosDataMap; }

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

    private EnsemblGeneData getGeneData(final String gene, boolean byName)
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

    public void populateGeneIdList(final List<String> uniqueGeneIds, final String chromosome, int position, int upstreamDistance)
    {
        // find the unique set of geneIds
        final List<EnsemblGeneData> matchedGenes = findGeneRegions(chromosome, position, upstreamDistance);

        for (final EnsemblGeneData geneData : matchedGenes)
        {
            if(!uniqueGeneIds.contains(geneData.GeneId))
                uniqueGeneIds.add(geneData.GeneId);
        }
    }

    public List<GeneAnnotation> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, int position,
            byte orientation, int upstreamDistance)
    {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        final List<EnsemblGeneData> matchedGenes = findGeneRegions(chromosome, position, upstreamDistance);

        // now look up relevant transcript and exon information
        for(final EnsemblGeneData geneData : matchedGenes)
        {
            final List<TranscriptData> transcriptDataList = mTranscriptDataMap.get(geneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            GeneAnnotation currentGene = new GeneAnnotation(svId, isStart, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.KaryotypeBand);

            currentGene.setGeneData(geneData);
            currentGene.setPositionalData(chromosome, position, orientation);

            // collect up all the relevant exons for each unique transcript to analyse as a collection
            for(TranscriptData transData : transcriptDataList)
            {
                Transcript transcript = extractTranscriptExonData(transData, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);

                    setAlternativeTranscriptPhasings(transcript, transData.exons(), position, orientation);

                    // annotate with preceding gene info if the up distance isn't set
                    if(!transcript.hasPrevSpliceAcceptorDistance())
                    {
                        setPrecedingGeneDistance(transcript, position);
                    }
                }
            }

            if(currentGene.transcripts().isEmpty() && mDownstreamGeneAnnotations.containsKey(geneData))
            {
                final TranscriptData trans = transcriptDataList.stream().filter(x -> x.IsCanonical).findAny().orElse(null);

                if(trans != null)
                {
                    int totalCodingBases = 0;
                    int codingBases = 0;
                    if(trans.CodingStart != null && trans.CodingEnd != null)
                    {
                        int[] codingData = Transcript.calcCodingBases(trans.CodingStart, trans.CodingEnd, trans.exons(), position);
                        totalCodingBases = codingData[TOTAL_CODING_BASES];
                    }

                    final Transcript postGeneTrans = new Transcript(
                            currentGene, trans.TransId, trans.TransName, trans.exons().size(), -1, trans.exons().size(),
                            -1, codingBases, totalCodingBases, trans.exons().size(), true,
                            trans.TransStart, trans.TransEnd, trans.CodingStart, trans.CodingEnd);

                    postGeneTrans.setBioType(trans.BioType);

                    currentGene.addTranscript(postGeneTrans);
                }
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    public List<GeneAnnotation> findGeneAnnotationsByOverlap(int svId, final String chromosome, int posStart, int posEnd)
    {
        // create gene and transcript data for any gene fully overlapped by the SV
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        final List<EnsemblGeneData> chrGeneList = mChrGeneDataMap.get(chromosome);

        if(chrGeneList == null)
            return geneAnnotations;

        for(final EnsemblGeneData geneData : chrGeneList)
        {
            if(!(posStart < geneData.GeneStart && posEnd > geneData.GeneEnd))
                continue;

            GeneAnnotation currentGene = new GeneAnnotation(svId, true, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.KaryotypeBand);

            currentGene.setGeneData(geneData);

            final TranscriptData transcriptData = mTranscriptDataMap.get(geneData.GeneId).stream()
                    .filter(x -> x.IsCanonical)
                    .findFirst().orElse(null);

            if (transcriptData == null)
                continue;

            Transcript transcript = extractTranscriptExonData(transcriptData, transcriptData.TransStart, currentGene);

            if(transcript != null)
            {
                currentGene.addTranscript(transcript);
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    private void setPrecedingGeneDistance(Transcript transcript, int position)
    {
        // annotate with preceding gene info if the up distance isn't set
        int precedingGeneSAPos = findPrecedingGeneSpliceAcceptorPosition(transcript.TransId);

        if(precedingGeneSAPos >= 0)
        {
            // if the breakend is after (higher for +ve strand) the nearest preceding splice acceptor, then the distance will be positive
            // and mean that the transcript isn't interrupted when used in a downstream fusion
            int preDistance = transcript.gene().Strand == 1 ? position - precedingGeneSAPos : precedingGeneSAPos - position;
            transcript.setSpliceAcceptorDistance(true, preDistance);
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

    public final List<EnsemblGeneData> findGenesByRegion(final String chromosome, int posStart, int posEnd)
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

    private List<EnsemblGeneData> findGeneRegions(final String chromosome, int position, int upstreamDistance)
    {
        final List<EnsemblGeneData> matchedGenes = Lists.newArrayList();

        final List<EnsemblGeneData> geneDataList = mChrGeneDataMap.get(chromosome);

        if(geneDataList == null)
            return matchedGenes;

        for(final EnsemblGeneData geneData : geneDataList)
        {
            int geneStartRange = geneData.Strand == 1 ? geneData.GeneStart - upstreamDistance : geneData.GeneStart;
            int geneEndRange = geneData.Strand == 1 ? geneData.GeneEnd : geneData.GeneEnd + upstreamDistance;

            if(position >= geneStartRange && position <= geneEndRange)
            {
                matchedGenes.add(geneData);
            }
        }

        for(Map.Entry<EnsemblGeneData,Integer> entry : mDownstreamGeneAnnotations.entrySet())
        {
            final EnsemblGeneData geneData = entry.getKey();

            if(matchedGenes.contains(geneData) || !geneData.Chromosome.equals(chromosome))
                continue;

            if((geneData.Strand == POS_STRAND && position >= geneData.GeneEnd && position <= geneData.GeneEnd + entry.getValue())
            || (geneData.Strand == NEG_STRAND && position <= geneData.GeneStart && position >= geneData.GeneStart - entry.getValue()))
            {
                matchedGenes.add(geneData);
            }
        }

        return matchedGenes;
    }

    public int findPrecedingGeneSpliceAcceptorPosition(int transId)
    {
        if(mTransSpliceAcceptorPosDataMap.isEmpty())
            return -1;

        Integer spliceAcceptorPos = mTransSpliceAcceptorPosDataMap.get(transId);
        return spliceAcceptorPos != null ? spliceAcceptorPos : -1;
    }

    public static Transcript extractTranscriptExonData(
            final TranscriptData transData, int position, final GeneAnnotation geneAnnotation)
    {
        final List<ExonData> exonList = transData.exons();

        if(exonList.isEmpty())
            return null;

        int exonMax = exonList.size();

        boolean isForwardStrand = geneAnnotation.Strand == 1;

        int upExonRank = -1;
        int upExonPhase = -1;
        int downExonRank = -1;
        int downExonPhase = -1;
        int nextUpDistance = -1;
        int nextDownDistance = -1;
        boolean isCodingTypeOverride = false;

        // first check for a position outside the exon boundaries
        final ExonData firstExon = exonList.get(0);
        final ExonData lastExon = exonList.get(exonList.size()-1);

        // for forward-strand transcripts the current exon is downstream, the previous is upstream
        // and the end-phase is taken from the upstream previous exon, the phase from the current downstream exon

        // for reverse-strand transcripts the current exon is upstream, the previous is downstream
        // and the end-phase is taken from the upstream (current) exon, the phase from the downstream (previous) exon

        // for each exon, the 'phase' is always the phase at the start of the exon in the direction of transcription
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
                        // the current exon is downstream, the previous one is upstream
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
                        // the previous exon in the list has the higher rank and is downstream
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

        int codingBases = 0;
        int totalCodingBases = 0;

        if(transData.CodingStart != null && transData.CodingEnd != null)
        {
            int[] codingData = Transcript.calcCodingBases(transData.CodingStart, transData.CodingEnd, exonList, position);
            totalCodingBases = codingData[TOTAL_CODING_BASES];
            codingBases = isForwardStrand ? codingData[CODING_BASES] : totalCodingBases - codingData[CODING_BASES];
        }

        Transcript transcript = new Transcript(geneAnnotation, transData.TransId, transData.TransName,
                upExonRank, upExonPhase, downExonRank, downExonPhase,
                codingBases, totalCodingBases,
                exonMax, transData.IsCanonical, transData.TransStart, transData.TransEnd,
                transData.CodingStart != null ? transData.CodingStart : null,
                transData.CodingEnd != null ? transData.CodingEnd : null);

        transcript.setBioType(transData.BioType);

        // if not set, leave the previous exon null and it will be taken from the closest upstream gene
        transcript.setSpliceAcceptorDistance(true, nextUpDistance >= 0 ? nextUpDistance : null);
        transcript.setSpliceAcceptorDistance(false, nextDownDistance >= 0 ? nextDownDistance : null);

        if(isCodingTypeOverride)
            transcript.setCodingType(CODING);

        return transcript;
    }

    public static int EXON_RANK_MIN = 0;
    public static int EXON_RANK_MAX = 1;
    public static int EXON_PHASE_MIN = 2;
    public static int EXON_PHASE_MAX = 3;

    public int[] getExonRankings(final String geneId, int position)
    {
        // finds the exon before and after this position, setting to -1 if before the first or beyond the last exon
        int[] exonData = new int[EXON_PHASE_MAX + 1];

        final TranscriptData transData = getTranscriptData(geneId, "");

        if (transData == null || transData.exons().isEmpty())
            return exonData;

        return getExonRankings(transData.Strand, transData.exons(), position);
    }

    public static int[] getExonRankings(int strand, final List<ExonData> exonDataList, int position)
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
            int position, byte orientation)
    {
        // collect exon phasings before the position on the upstream and after it on the downstream
        boolean isUpstream = (transcript.gene().Strand * orientation) > 0;
        boolean forwardStrand = (transcript.gene().Strand == 1);

        Map<Integer,Integer> alternativePhasing = Maps.newHashMap();

        int transPhase = isUpstream ? transcript.ExonUpstreamPhase : transcript.ExonDownstreamPhase;
        int transRank = isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream;

        for (ExonData exonData : exonDataList)
        {

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
            int exonsSkipped;

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

    public boolean load(boolean delayTranscriptLoading)
    {
        if(!loadEnsemblGeneData(mDataPath, mRestrictedGeneIdList, mChrGeneDataMap, mRefGenomeVersion))
            return false;

        if(!delayTranscriptLoading)
        {
            if(!EnsemblDataLoader.loadTranscriptData(mDataPath, mTranscriptDataMap, mRestrictedGeneIdList, mRequireExons, mCanonicalTranscriptsOnly))
                return false;

            if(mRequireProteinDomains && !loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, Lists.newArrayList()))
                return false;

            if(mRequireSplicePositions && !loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, Lists.newArrayList()))
                return false;
        }

        return true;
    }

    public boolean loadTranscriptData(final List<String> restrictedGeneIds)
    {
        if(!EnsemblDataLoader.loadTranscriptData(mDataPath, mTranscriptDataMap, restrictedGeneIds, mRequireExons, mCanonicalTranscriptsOnly))
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

        if(mRequireProteinDomains && !loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, uniqueTransIds))
            return false;

        if(mRequireSplicePositions && !loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, uniqueTransIds))
            return false;

        return true;
    }

    public static Integer[] getProteinDomainPositions(final TranscriptProteinData proteinData, final TranscriptData transData)
    {
        Integer[] domainPositions = {null, null};

        if(transData.exons().isEmpty())
            return domainPositions;

        Integer codingStart = transData.CodingStart;
        Integer codingEnd = transData.CodingEnd;

        if(codingStart == null || codingEnd == null)
            return domainPositions;

        int preProteinBases = proteinData.SeqStart * 3;
        int proteinBases = (proteinData.SeqEnd - proteinData.SeqStart) * 3;

        int proteinStart = -1;
        int proteinEnd = -1;

        if(transData.Strand == 1)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exonData = transData.exons().get(i);

                if(exonData.ExonEnd < codingStart)
                    continue;

                if(preProteinBases > 0)
                {
                    int refStartPos = max(codingStart, exonData.ExonStart);
                    int exonCodingBases = exonData.ExonEnd - refStartPos;

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

                int startPos = max(exonData.ExonStart, proteinStart);
                int exonBases = exonData.ExonEnd - startPos;

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
                    int refStartPos = min(codingEnd, exonData.ExonEnd);
                    int exonCodingBases = refStartPos - exonData.ExonStart;

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

                int startPos = min(exonData.ExonEnd, proteinEnd);
                int exonBases = startPos - exonData.ExonStart;

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
