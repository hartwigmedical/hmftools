package com.hartwig.hmftools.common.ensemblcache;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptProteinData;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.loadTranscriptSpliceAcceptorData;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.gene.TranscriptProteinData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class EnsemblDataCache
{
    private final String mDataPath;
    private final RefGenomeVersion mRefGenomeVersion;

    private final Map<String,List<TranscriptData>> mTranscriptByGeneIdMap; // transcripts keyed by geneId
    private final Map<Integer,TranscriptData> mTranscriptsByTransIdMap; // transcripts keyed by transId
    private final Map<String,List<GeneData>> mChrGeneDataMap; // genes keyed by chromosome
    private final Map<Integer,List<TranscriptProteinData>> mEnsemblProteinDataMap;
    private final Map<Integer,Integer> mTransSpliceAcceptorPosDataMap;
    private final Map<String,GeneData> mGeneDataMap; // keyed by geneId
    private final Map<String,GeneData> mGeneNameIdMap; // for faster look-up by name

    // whether to load more details information for each transcript - exons, protein domains, splice positions etc
    private boolean mRequireExons;
    private boolean mRequireProteinDomains;
    private boolean mRequireSplicePositions;
    private boolean mCanonicalTranscriptsOnly;
    private boolean mRequireGeneSynonyms;
    private boolean mRequireNonEnsemblTranscripts;

    private final Map<GeneData,Integer> mDownstreamGeneAnnotations;
    private final List<GeneData> mAlternativeGeneData;
    private final List<String> mRestrictedGeneIdList = Lists.newArrayList();

    public static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";
    public static final String ENSEMBL_DATA_DIR_CFG = "Ensembl data file directory";

    public EnsemblDataCache(final ConfigBuilder configBuilder)
    {
        this(configBuilder.getValue(ENSEMBL_DATA_DIR), RefGenomeVersion.from(configBuilder));
    }

    public EnsemblDataCache(final String dataPath, final RefGenomeVersion refGenomeVersion)
    {
        mDataPath = checkAddDirSeparator(dataPath);
        mRefGenomeVersion = refGenomeVersion;

        mTranscriptByGeneIdMap = Maps.newHashMap();
        mTranscriptsByTransIdMap = Maps.newHashMap();
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
        mRequireNonEnsemblTranscripts = false;
        mDownstreamGeneAnnotations = Maps.newHashMap();
        mAlternativeGeneData = Lists.newArrayList();
    }

    public static void addEnsemblDir(final ConfigBuilder configBuilder)
    {
        addEnsemblDir(configBuilder, false);
    }

    public static void addEnsemblDir(final ConfigBuilder configBuilder, boolean required)
    {
        configBuilder.addPath(ENSEMBL_DATA_DIR, required, ENSEMBL_DATA_DIR_CFG);
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
    public boolean hasDownstreamGeneAnnotation(final GeneData geneData) { return mDownstreamGeneAnnotations.containsKey(geneData); }

    public List<GeneData> getAlternativeGeneData() { return mAlternativeGeneData; }

    public void setRequiredData(boolean exons, boolean proteinDomains, boolean splicePositions, boolean canonicalOnly)
    {
        mRequireExons = exons;
        mRequireSplicePositions = splicePositions;
        mRequireProteinDomains = proteinDomains;
        mCanonicalTranscriptsOnly = canonicalOnly;
    }

    public void setRequireGeneSynonyms() { mRequireGeneSynonyms = true; }
    public void setRequireNonEnsemblTranscripts() { mRequireNonEnsemblTranscripts = true; }

    public Map<String,List<TranscriptData>> getTranscriptDataMap() { return mTranscriptByGeneIdMap; }
    public Map<String,List<GeneData>> getChrGeneDataMap() { return mChrGeneDataMap; }
    public Map<Integer,List<TranscriptProteinData>> getTranscriptProteinDataMap() { return mEnsemblProteinDataMap; }

    public final GeneData getGeneDataByName(final String geneName)
    {
        if(!mGeneNameIdMap.isEmpty())
            return mGeneNameIdMap.get(geneName);

        return getGeneData(geneName, true);
    }

    public GeneData getGeneDataById(final String geneId)
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

    public void createTranscriptIdMap()
    {
        if(!mTranscriptsByTransIdMap.isEmpty())
            return;

        for(List<TranscriptData> transcriptDataList : mTranscriptByGeneIdMap.values())
        {
            transcriptDataList.forEach(x -> mTranscriptsByTransIdMap.put(x.TransId, x));
        }
    }

    public List<TranscriptData> getTranscripts(final String geneId)
    {
        return mTranscriptByGeneIdMap.get(geneId);
    }

    public TranscriptData getCanonicalTranscriptData(final String geneId) { return getTranscriptData(geneId, ""); }

    public TranscriptData getTranscriptData(final String geneId, final String transcriptId)
    {
        // leave transcriptId empty to retrieve the canonical transcript
        final List<TranscriptData> transDataList = mTranscriptByGeneIdMap.get(geneId);

        if(transDataList == null || transDataList.isEmpty())
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

    public TranscriptData getTranscriptData(final int transId) { return mTranscriptsByTransIdMap.get(transId); }

    public final List<GeneData> findGenesByRegion(final String chromosome, int posStart, int posEnd)
    {
        // find genes if any of their transcripts are within this position
        List<GeneData> genesList = Lists.newArrayList();

        final List<GeneData> geneDataList = mChrGeneDataMap.get(chromosome);

        for(final GeneData geneData : geneDataList)
        {
            if(posStart > geneData.GeneEnd || posEnd < geneData.GeneStart)
                continue;

            final List<TranscriptData> transList = mTranscriptByGeneIdMap.get(geneData.GeneId);

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

    public List<GeneData> findGeneRegions(final String chromosome, int position, int upstreamDistance)
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

    public boolean load(boolean delayTranscriptLoading)
    {
        if(!loadEnsemblGeneData(mDataPath, mRestrictedGeneIdList, mChrGeneDataMap, mRefGenomeVersion, mRequireGeneSynonyms))
            return false;

        if(!delayTranscriptLoading)
        {
            if(!EnsemblDataLoader.loadTranscriptData(
                    mDataPath, mTranscriptByGeneIdMap, mRestrictedGeneIdList, mRequireExons, mCanonicalTranscriptsOnly,
                    mRequireNonEnsemblTranscripts, Collections.emptyList()))
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
        return loadTranscriptData(restrictedGeneIds, Collections.emptyList());
    }

    public boolean loadTranscriptData(final List<String> restrictedGeneIds, final List<String> nonCanonicalTrans)
    {
        if(!EnsemblDataLoader.loadTranscriptData(
                mDataPath, mTranscriptByGeneIdMap, restrictedGeneIds, mRequireExons, mCanonicalTranscriptsOnly,
                mRequireNonEnsemblTranscripts, nonCanonicalTrans))
        {
            return false;
        }

        Set<Integer> uniqueTransIds = Sets.newHashSet();

        for(List<TranscriptData> transDataList : mTranscriptByGeneIdMap.values())
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
