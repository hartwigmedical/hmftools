package com.hartwig.hmftools.common.ensemblcache;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptProteinData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptSpliceAcceptorData;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.CodingBaseData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptProteinData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class EnsemblDataCache
{
    private final String mDataPath;
    private final RefGenomeVersion mRefGenomeVersion;

    private final Map<String,List<TranscriptData>> mTranscriptDataMap; // transcripts keyed by geneId
    private final Map<String,List<GeneData>> mChrGeneDataMap; // genes keyed by chromosome
    private final Map<Integer,List<TranscriptProteinData>> mEnsemblProteinDataMap;
    private final Map<Integer,Integer> mTransSpliceAcceptorPosDataMap;
    private final Map<String, GeneData> mGeneDataMap; // keyed by geneId
    private final Map<String, GeneData> mGeneNameIdMap; // for faster look-up by name

    private GeneNameMapping mGeneNameMapping;

    // whether to load more details information for each transcript - exons, protein domains, splice positions etc
    private boolean mRequireExons;
    private boolean mRequireProteinDomains;
    private boolean mRequireSplicePositions;
    private boolean mCanonicalTranscriptsOnly;
    private boolean mRequireGeneSynonyms;

    private final Map<GeneData,Integer> mDownstreamGeneAnnotations;
    private final List<GeneData> mAlternativeGeneData;
    private final List<String> mRestrictedGeneIdList = Lists.newArrayList();

    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir"; // eventually deprecated
    public static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";
    public static final String ENSEMBL_DATA_DIR_CFG = "Ensembl data file directory";

    public EnsemblDataCache(final CommandLine cmd, final RefGenomeVersion refGenomeVersion)
    {
        this(cmd.hasOption(ENSEMBL_DATA_DIR) ?
                cmd.getOptionValue(ENSEMBL_DATA_DIR) : cmd.getOptionValue(GENE_TRANSCRIPTS_DIR), refGenomeVersion);
    }

    public EnsemblDataCache(final String dataPath, final RefGenomeVersion refGenomeVersion)
    {
        mDataPath = checkAddDirSeparator(dataPath);
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
        mRequireGeneSynonyms = false;
        mDownstreamGeneAnnotations = Maps.newHashMap();
        mAlternativeGeneData = Lists.newArrayList();
        mGeneNameMapping = null;
    }

    public static void addEnsemblDir(final Options options)
    {
        options.addOption(ENSEMBL_DATA_DIR, true, ENSEMBL_DATA_DIR_CFG);
        options.addOption(GENE_TRANSCRIPTS_DIR, true, ENSEMBL_DATA_DIR_CFG);
    }

    public void setRestrictedGeneIdList(final List<String> geneIds)
    {
        mRestrictedGeneIdList.clear();
        mRestrictedGeneIdList.addAll(geneIds);
    }

    public void addDownstreamGeneAnnotations(final GeneData geneData, int distance)
    {
        mDownstreamGeneAnnotations.put(geneData, distance);
    }

    public final List<GeneData> getAlternativeGeneData() { return mAlternativeGeneData; }

    public void setRequiredData(boolean exons, boolean proteinDomains, boolean splicePositions, boolean canonicalOnly)
    {
        mRequireExons = exons;
        mRequireSplicePositions = splicePositions;
        mRequireProteinDomains = proteinDomains;
        mCanonicalTranscriptsOnly = canonicalOnly;
    }

    public void setRequireGeneSynonyms() { mRequireGeneSynonyms = true; }

    public final Map<String,List<TranscriptData>> getTranscriptDataMap() { return mTranscriptDataMap; }
    public final Map<String,List<GeneData>> getChrGeneDataMap() { return mChrGeneDataMap; }
    public Map<Integer,List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }
    public Map<Integer,Integer> getTransSpliceAcceptorPosDataMap() { return mTransSpliceAcceptorPosDataMap; }

    public final GeneData getGeneDataByName(final String geneName)
    {
        if(!mGeneNameIdMap.isEmpty())
            return mGeneNameIdMap.get(geneName);

        return getGeneData(geneName, true);
    }

    public final GeneData getGeneDataById(final String geneId)
    {
        if(!mGeneDataMap.isEmpty())
            return mGeneDataMap.get(geneId);

        return getGeneData(geneId, false);
    }

    private GeneData getGeneData(final String gene, boolean byName)
    {
        for(Map.Entry<String, List<GeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final GeneData geneData : entry.getValue())
            {
                if((byName && geneData.GeneName.equals(gene)) || (!byName && geneData.GeneId.equals(gene)))
                    return geneData;
            }
        }

        return null;
    }

    public GeneData getGeneDataBySynonym(final String geneName, final String chromosome)
    {
        final List<GeneData> geneDataList = mChrGeneDataMap.get(chromosome);
        if(geneDataList == null)
            return null;

        return geneDataList.stream().filter(x -> x.hasSynonym(geneName)).findFirst().orElse(null);
    }

    public void createGeneIdDataMap()
    {
        if(!mGeneDataMap.isEmpty())
            return;

        for(Map.Entry<String, List<GeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final GeneData geneData : entry.getValue())
            {
                mGeneDataMap.put(geneData.GeneId, geneData);
            }
        }
    }

    public void createGeneNameIdMap()
    {
        if(!mGeneNameIdMap.isEmpty())
            return;

        for(Map.Entry<String, List<GeneData>> entry : mChrGeneDataMap.entrySet())
        {
            for(final GeneData geneData : entry.getValue())
            {
                mGeneNameIdMap.put(geneData.GeneName, geneData);
            }
        }
    }

    public GeneNameMapping getGeneMappings()
    {
        if(mGeneNameMapping == null)
            mGeneNameMapping = new GeneNameMapping();

        return mGeneNameMapping;
    }

    public List<TranscriptData> getTranscripts(final String geneId)
    {
        return mTranscriptDataMap.get(geneId);
    }

    public void populateGeneIdList(final List<String> uniqueGeneIds, final String chromosome, int position, int upstreamDistance)
    {
        // find the unique set of geneIds
        final List<GeneData> matchedGenes = findGeneRegions(chromosome, position, upstreamDistance);

        for (final GeneData geneData : matchedGenes)
        {
            if(!uniqueGeneIds.contains(geneData.GeneId))
                uniqueGeneIds.add(geneData.GeneId);
        }
    }

    public List<BreakendGeneData> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, int position,
            byte orientation, int upstreamDistance)
    {
        List<BreakendGeneData> geneAnnotations = Lists.newArrayList();

        final List<GeneData> matchedGenes = findGeneRegions(chromosome, position, upstreamDistance);

        // now look up relevant transcript and exon information
        for(final GeneData geneData : matchedGenes)
        {
            final List<TranscriptData> transcriptDataList = mTranscriptDataMap.get(geneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            BreakendGeneData currentGene = new BreakendGeneData(svId, isStart, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.KaryotypeBand);

            currentGene.setGeneData(geneData);
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
                        setPrecedingGeneDistance(transcript, position);
                    }
                }
            }

            if(currentGene.transcripts().isEmpty() && mDownstreamGeneAnnotations.containsKey(geneData))
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

        final List<GeneData> altMappingGenes = mAlternativeGeneData.stream()
                .filter(x -> x.Chromosome.equals(chromosome))
                .filter(x -> positionWithin(position, x.GeneStart, x.GeneEnd))
                .collect(Collectors.toList());

        for(final GeneData altGeneData : altMappingGenes)
        {
            final List<TranscriptData> transcriptDataList = mTranscriptDataMap.get(altGeneData.GeneId);

            if (transcriptDataList == null || transcriptDataList.isEmpty())
                continue;

            final TranscriptData trans = transcriptDataList.stream().filter(x -> x.IsCanonical).findFirst().orElse(null);

            BreakendGeneData geneAnnotation = new BreakendGeneData(svId, isStart, altGeneData.GeneName, altGeneData.GeneId,
                    altGeneData.Strand, altGeneData.KaryotypeBand);

            geneAnnotation.setGeneData(altGeneData);

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

    private void setPrecedingGeneDistance(BreakendTransData transcript, int position)
    {
        // annotate with preceding gene info if the up distance isn't set
        int precedingGeneSAPos = findPrecedingGeneSpliceAcceptorPosition(transcript.transId());

        if(precedingGeneSAPos >= 0)
        {
            // if the breakend is after (higher for +ve strand) the nearest preceding splice acceptor, then the distance will be positive
            // and mean that the transcript isn't interrupted when used in a downstream fusion
            int preDistance = transcript.gene().Strand == POS_STRAND ? position - precedingGeneSAPos : precedingGeneSAPos - position;
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

    public final List<GeneData> findGenesByRegion(final String chromosome, int posStart, int posEnd)
    {
        // find genes if any of their transcripts are within this position
        List<GeneData> genesList = Lists.newArrayList();

        final List<GeneData> geneDataList = mChrGeneDataMap.get(chromosome);

        for(final GeneData geneData : geneDataList)
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

    private List<GeneData> findGeneRegions(final String chromosome, int position, int upstreamDistance)
    {
        final List<GeneData> matchedGenes = Lists.newArrayList();

        final List<GeneData> geneDataList = mChrGeneDataMap.get(chromosome);

        if(geneDataList == null)
            return matchedGenes;

        for(final GeneData geneData : geneDataList)
        {
            int geneStartRange = geneData.Strand == 1 ? geneData.GeneStart - upstreamDistance : geneData.GeneStart;
            int geneEndRange = geneData.Strand == 1 ? geneData.GeneEnd : geneData.GeneEnd + upstreamDistance;

            if(position >= geneStartRange && position <= geneEndRange)
            {
                matchedGenes.add(geneData);
            }
        }

        for(Map.Entry<GeneData,Integer> entry : mDownstreamGeneAnnotations.entrySet())
        {
            final GeneData geneData = entry.getKey();

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

    public static BreakendTransData createBreakendTranscriptData(
            final TranscriptData transData, int position, final BreakendGeneData geneAnnotation)
    {
        final List<ExonData> exonList = transData.exons();

        if(exonList.isEmpty())
            return null;

        boolean isForwardStrand = geneAnnotation.Strand == POS_STRAND;
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

    public static final int EXON_RANK_MIN = 0;
    public static final int EXON_RANK_MAX = 1;

    public int[] getExonRankings(final String geneId, int position)
    {
        // finds the exon before and after this position, setting to -1 if before the first or beyond the last exon
        int[] exonData = new int[EXON_RANK_MAX + 1];

        final TranscriptData transData = getTranscriptData(geneId, "");

        if (transData == null || transData.exons().isEmpty())
            return exonData;

        return getExonRankings(transData.Strand, transData.exons(), position);
    }

    public static int[] getExonRankings(int strand, final List<ExonData> exonDataList, int position)
    {
        int[] exonData = new int[EXON_RANK_MAX + 1];

        // first test a position outside the range of the exons
        final ExonData firstExon = exonDataList.get(0);
        final ExonData lastExon = exonDataList.get(exonDataList.size() - 1);

        if((position < firstExon.Start && strand == POS_STRAND) || (position > lastExon.End && strand == NEG_STRAND))
        {
            // before the start of the transcript
            exonData[EXON_RANK_MIN] = 0;
            exonData[EXON_RANK_MAX] = 1;
        }
        else if((position < firstExon.Start && strand == NEG_STRAND) || (position > lastExon.End && strand == POS_STRAND))
        {
            // past the end of the transcript
            exonData[EXON_RANK_MIN] = exonDataList.size();
            exonData[EXON_RANK_MAX] = -1;
        }
        else
        {
            for(int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData transExonData = exonDataList.get(i);
                final ExonData nextTransExonData = i < exonDataList.size() - 1 ? exonDataList.get(i+1) : null;

                if(position == transExonData.End || position == transExonData.Start)
                {
                    // position matches the bounds of an exon
                    exonData[EXON_RANK_MIN] = transExonData.Rank;
                    exonData[EXON_RANK_MAX] = transExonData.Rank;
                    break;
                }

                if(position >= transExonData.Start && position <= transExonData.End)
                {
                    // position matches within or at the bounds of an exon
                    exonData[EXON_RANK_MIN] = transExonData.Rank;
                    exonData[EXON_RANK_MAX] = transExonData.Rank;
                    break;
                }

                if(nextTransExonData != null && position > transExonData.End && position < nextTransExonData.Start)
                {
                    if(strand == 1)
                    {
                        exonData[EXON_RANK_MIN] = transExonData.Rank;
                        exonData[EXON_RANK_MAX] = nextTransExonData.Rank;
                    }
                    else
                    {
                        exonData[EXON_RANK_MIN] = nextTransExonData.Rank;
                        exonData[EXON_RANK_MAX] = transExonData.Rank;
                    }

                    break;
                }
            }
        }

        return exonData;
    }

    public static void setAlternativeTranscriptPhasings(BreakendTransData transcript, final List<ExonData> exonDataList, int position, byte orientation)
    {
        // collect exon phasings before the position on the upstream and after it on the downstream
        boolean isUpstream = (transcript.gene().Strand * orientation) > 0;
        boolean forwardStrand = (transcript.gene().Strand == POS_STRAND);

        Map<Integer,Integer> alternativePhasing = transcript.getAlternativePhasing();
        alternativePhasing.clear();

        int transPhase = isUpstream ? transcript.Phase : transcript.Phase;
        int transRank = isUpstream ? transcript.ExonUpstream : transcript.ExonDownstream;

        for (ExonData exonData : exonDataList)
        {
            if(isUpstream == forwardStrand)
            {
                if (exonData.Start > position || transRank == exonData.Rank)
                    break;
            }
            else
            {
                if (position > exonData.End || transRank == exonData.Rank)
                    continue;
            }

            int exonPhase = isUpstream ? exonData.PhaseEnd : exonData.PhaseStart;
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

    public boolean load(boolean delayTranscriptLoading)
    {
        if(!loadEnsemblGeneData(mDataPath, mRestrictedGeneIdList, mChrGeneDataMap, mRefGenomeVersion, mRequireGeneSynonyms))
            return false;

        if(!delayTranscriptLoading)
        {
            if(!EnsemblDataLoader.loadTranscriptData(
                    mDataPath, mTranscriptDataMap, mRestrictedGeneIdList, mRequireExons, mCanonicalTranscriptsOnly, Lists.newArrayList()))
            {
                return false;
            }

            if(mRequireProteinDomains && !loadTranscriptProteinData(mDataPath, mEnsemblProteinDataMap, Sets.newHashSet()))
                return false;

            if(mRequireSplicePositions && !loadTranscriptSpliceAcceptorData(mDataPath, mTransSpliceAcceptorPosDataMap, Sets.newHashSet()))
                return false;
        }

        return true;
    }

    public boolean loadTranscriptData(final List<String> restrictedGeneIds)
    {
        return loadTranscriptData(restrictedGeneIds, Lists.newArrayList());
    }

    public boolean loadTranscriptData(final List<String> restrictedGeneIds, final List<String> nonCanonicalTrans)
    {
        if(!EnsemblDataLoader.loadTranscriptData(
                mDataPath, mTranscriptDataMap, restrictedGeneIds, mRequireExons, mCanonicalTranscriptsOnly, nonCanonicalTrans))
        {
            return false;
        }

        Set<Integer> uniqueTransIds = Sets.newHashSet();

        for(List<TranscriptData> transDataList : mTranscriptDataMap.values())
        {
            transDataList.forEach(x -> uniqueTransIds.add(x.TransId));
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

                if(exonData.End < codingStart)
                    continue;

                if(preProteinBases > 0)
                {
                    int refStartPos = max(codingStart, exonData.Start);
                    int exonCodingBases = exonData.End - refStartPos;

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

                int startPos = max(exonData.Start, proteinStart);
                int exonBases = exonData.End - startPos;

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

                if(exonData.Start > codingEnd)
                    continue;

                if(preProteinBases > 0)
                {
                    int refStartPos = min(codingEnd, exonData.End);
                    int exonCodingBases = refStartPos - exonData.Start;

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

                int startPos = min(exonData.End, proteinEnd);
                int exonBases = startPos - exonData.Start;

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
