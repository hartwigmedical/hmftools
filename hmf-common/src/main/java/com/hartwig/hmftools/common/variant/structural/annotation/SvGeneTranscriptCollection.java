package com.hartwig.hmftools.common.variant.structural.annotation;

import static java.lang.Math.abs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<Integer, List<GeneAnnotation>> mSvIdGeneTranscriptsMap;
    private SortedSetMultimap<String, HmfTranscriptRegion> mGenesByChromosomeMap;
    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private Map<String, List<TranscriptExonData>> mGeneTransExonDataMap;
    private Map<String, EnsemblGeneData> mEnsemblGeneDataMap;

    // to get a wider range of candidate genes and filter by promotor distance later on
    public static int PRE_GENE_PROMOTOR_DISTANCE = 100000;

    public static String SV_GENE_TRANSCRIPTS_FILE_SUFFIX = "sv_ensembl_data.csv";

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mSvIdGeneTranscriptsMap = new HashMap();
        mGeneTransExonDataMap = new HashMap();
        mEnsemblGeneDataMap = new HashMap();
        mGenesByChromosomeMap = null;
        mAllGenesMap = null;
    }

    public final Map<Integer, List<GeneAnnotation>> getSvIdGeneTranscriptsMap() { return mSvIdGeneTranscriptsMap; }
    public final Map<String, List<TranscriptExonData>> getGeneExonDataMap() { return mGeneTransExonDataMap; }

    public void setDataPath(final String dataPath)
    {
        mDataPath = dataPath;
    }

    // Gene,CanonicalTranscriptId,Strand,TransId,Trans,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,
    // ExonPhase,ExonEndPhase,CodingStart,CodingEnd
    private static int TE_GENE_ID = 0;
    private static int TE_CANONICAL = 1;
    private static int TE_STRAND = 2;
    private static int TE_TRANS_ID = 3;
    private static int TE_TRANS_NAME = 4;
    private static int TE_BIOTYPE = 5;
    private static int TE_TRANS_START = 6;
    private static int TE_TRANS_END = 7;
    private static int TE_EXON_RANK = 8;
    private static int TE_EXON_START = 9;
    private static int TE_EXON_END = 10;
    private static int TE_PHASE = 11;
    private static int TE_PHASE_END = 12;
    private static int TE_CODING_START = 13;
    private static int TE_CODING_END = 14;

    public boolean loadTranscriptExonData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl gene-exon data file({})", filename);
                return false;
            }

            int exonCount = 0;
            String currentGene = "";
            List<TranscriptExonData> transExonDataList = null;

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String geneId = items[TE_GENE_ID];

                if(!geneId.equals(currentGene))
                {
                    currentGene = geneId;
                    transExonDataList = Lists.newArrayList();
                    mGeneTransExonDataMap.put(geneId, transExonDataList);
                }

                // Gene,CanonicalTranscriptId,Strand,TransId,Trans,TransStart,TransEnd,ExonRank,ExonStart,ExonEnd,
                // ExonPhase,ExonEndPhase,CodingStart,CodingEnd

                Long codingStart = !items[TE_CODING_START].equals("NULL") ? Long.parseLong(items[TE_CODING_START]) : null;
                Long codingEnd = !items[TE_CODING_END].equals("NULL") ? Long.parseLong(items[TE_CODING_END]) : null;

                TranscriptExonData exonData = new TranscriptExonData(
                        geneId, items[TE_TRANS_NAME], Integer.parseInt(items[TE_TRANS_ID]),
                        Boolean.parseBoolean(items[TE_CANONICAL]), Byte.parseByte(items[TE_STRAND]),
                        Long.parseLong(items[TE_TRANS_START]), Long.parseLong(items[TE_TRANS_END]),
                        Long.parseLong(items[TE_EXON_START]), Long.parseLong(items[TE_EXON_END]),
                        Integer.parseInt(items[TE_EXON_RANK]), Integer.parseInt(items[TE_PHASE]), Integer.parseInt(items[TE_PHASE_END]),
                        codingStart, codingEnd, items[TE_BIOTYPE]);

                transExonDataList.add(exonData);
                ++exonCount;

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} gene records, {} exon", mGeneTransExonDataMap.size(), exonCount);

            mGenesByChromosomeMap = HmfGenePanelSupplier.allGenesPerChromosomeMap37();

            mAllGenesMap = Maps.newHashMap();
            for (final HmfTranscriptRegion region : mGenesByChromosomeMap.values())
            {
                mAllGenesMap.put(region.gene(), region);
            }

        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load sample gene annotations({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    // GeneId,GeneName,Chromosome,Strand,GeneStart,GeneEnd,EntrezIds,KaryotypeBand,Synonyms
    private static int GD_ID = 0;
    private static int GD_NAME = 1;
    private static int GD_CHR = 2;
    private static int GD_STRAND = 3;
    private static int GD_START = 4;
    private static int GD_END = 5;
    private static int GD_ENTREZ = 6;
    private static int GD_BAND = 7;
    private static int GD_SYN = 8;

    public boolean loadEnsemblGeneData(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty Ensembl gene data file({})", filename);
                return false;
            }

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                String[] items = line.split(",");

                final String geneId = items[GD_ID];

                EnsemblGeneData geneData = new EnsemblGeneData(
                        geneId, items[GD_NAME], items[GD_CHR], Byte.parseByte(items[GD_STRAND]),
                        Long.parseLong(items[GD_START]), Long.parseLong(items[GD_END]),
                        items[GD_ENTREZ], items[GD_BAND], items[GD_SYN]);

                mEnsemblGeneDataMap.put(geneId, geneData);

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} gene records", mEnsemblGeneDataMap.size());
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load Ensembl gene ({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public List<GeneAnnotation> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, long position, byte orientation)
    {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        final SortedSet<HmfTranscriptRegion> geneRegions = mGenesByChromosomeMap.get(chromosome);

        final List<HmfTranscriptRegion> matchedGenes = findGeneRegions(position, orientation, geneRegions);

        // now look up relevant transcript and exon information
        for(final HmfTranscriptRegion geneRegion : matchedGenes)
        {
            final List<TranscriptExonData> transExonDataList = mGeneTransExonDataMap.get(geneRegion.geneID());

            if (transExonDataList == null || transExonDataList.isEmpty())
                continue;

            final EnsemblGeneData geneData = mEnsemblGeneDataMap.get(geneRegion.geneID());

            if(geneData == null)
            {
                LOGGER.warn("gene({}) data not found in Ensembl cache", geneRegion.gene());
                continue;
            }

            GeneAnnotation currentGene = new GeneAnnotation(svId, isStart, geneRegion.gene(), geneData.GeneId,
                    geneData.Strand, geneData.Synonyms, geneData.EntrezIds, geneData.KaryotypeBand);

            for(int i = 0; i < transExonDataList.size(); ++i)
            {
                List<TranscriptExonData> transcriptExons = Lists.newArrayList();
                int transId = transExonDataList.get(i).TransId;

                int j = i;
                for(; j < transExonDataList.size(); ++j)
                {
                    if(transExonDataList.get(j).TransId != transId)
                        break;

                    transcriptExons.add(transExonDataList.get(j));
                }

                Transcript transcript = extractTranscript(transcriptExons, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);
                    // LOGGER.debug("gene({}) transcript({}) added info from exon records({} -> {})", geneData.GeneName, transcript.transcriptId(), i, j);
                }

                if(j == transExonDataList.size() - 1)
                    break;

                i = j - 1;
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    private Transcript extractTranscript(final List<TranscriptExonData> transcriptExons, long position, final GeneAnnotation geneAnnotation)
    {
        int exonMax = transcriptExons.size();

        final TranscriptExonData first = transcriptExons.get(0);

        boolean isCoding = first.CodingStart != null && first.CodingEnd != null;
        long codingStart = first.CodingStart != null ? first.CodingStart : 0;
        long codingEnd = first.CodingEnd != null ? first.CodingEnd : 0;
        boolean isForwardStrand = geneAnnotation.strand() == 1;

        // for the given position, determine how many coding bases occur prior to the position
        // in the direction of the transcript
        // strand direction will be corrected for afterwards

        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        long codingBases = 0;
        long totalCodingBases = 0;
        long transcriptStart = 0;
        long transcriptEnd = 0;

        // previous here will be the earlier exon, ordered by increasing position (ie regardless of strand direction)
        int prevExonRank = -1;
        int prevExonPhase = 0;
        int prevExonEndPhase = 0;

        // similarly the next exon will be exon immediately after the position for exons which increase with positino
        int nextExonRank = -1;
        int nextExonPhase = 0;
        int nextExonEndPhase = 0;

        for (int index = 0; index < transcriptExons.size(); ++index)
        {
            final TranscriptExonData exonData = transcriptExons.get(index);
            long exonStart = exonData.ExonStart;
            long exonEnd = exonData.ExonEnd;

            if(index == 0)
                transcriptStart = exonStart;

            if(index == transcriptExons.size() - 1)
                transcriptEnd = exonEnd;

            if(position >= exonStart && position <= exonEnd)
            {
                // falls within an exon
                prevExonRank = nextExonRank = exonData.ExonRank;

                if(isForwardStrand)
                {
                    prevExonEndPhase = exonData.ExonPhase;
                    nextExonPhase = exonData.ExonPhaseEnd;

                    // won't be used
                    prevExonPhase = prevExonEndPhase;
                    nextExonEndPhase = nextExonPhase;
                }
                else
                {
                    prevExonPhase = exonData.ExonPhaseEnd;
                    nextExonEndPhase = exonData.ExonPhase;

                    prevExonEndPhase = prevExonPhase;
                    nextExonPhase = nextExonEndPhase;
                }
            }
            else if(position > exonEnd)
            {
                // continue updating this until past the position
                prevExonRank = exonData.ExonRank;
                prevExonPhase = exonData.ExonPhase;
                prevExonEndPhase = exonData.ExonPhaseEnd;
            }
            else if(position < exonStart && nextExonRank == -1)
            {
                // set at the first exon past this position
                nextExonRank = exonData.ExonRank;
                nextExonPhase = exonData.ExonPhase;
                nextExonEndPhase = exonData.ExonPhaseEnd;
            }

            if(!isCoding)
                continue;

            if(!inCodingRegion)
            {
                if(exonEnd >= codingStart)
                {
                    // coding region begins in this exon
                    inCodingRegion = true;

                    totalCodingBases += exonEnd - codingStart + 1;

                    // check whether the position falls in this exon and if so before or after the coding start
                    if(position >= codingStart)
                    {
                        if(position < exonEnd)
                            codingBases += position - codingStart + 1;
                        else
                            codingBases += exonEnd - codingStart + 1;
                    }
                }
            }
            else if(!codingRegionEnded)
            {
                if(exonStart > codingEnd)
                {
                    codingRegionEnded = true;
                }
                else if(exonEnd > codingEnd)
                {
                    // coding region ends in this exon
                    codingRegionEnded = true;

                    totalCodingBases += codingEnd - exonStart + 1;

                    if(position >= exonStart)
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

                    if(position >= exonStart)
                    {
                        if (position < exonEnd)
                            codingBases += position - exonStart + 1;
                        else
                            codingBases += exonEnd - exonStart + 1;
                    }
                }
            }
        }

        if(prevExonRank == -1)
        {
            if(!isForwardStrand)
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                //                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on reverse strand",
                //                        transcriptStableId, position, nextExonRank, exonMax);
                return null;
            }
            else
            {
                prevExonRank = 0;
                prevExonPhase = -1;
                prevExonEndPhase = -1;
            }
        }
        else if(nextExonRank == -1)
        {
            if(isForwardStrand)
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                //                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on forward strand",
                //                        transcriptStableId, position, prevExonRank, exonMax);
                return null;
            }
            else
            {
                nextExonRank = 0;
                nextExonPhase = -1;
                nextExonEndPhase = -1;
            }
        }

        if(nextExonRank < 0 || prevExonRank < 0 || abs(nextExonRank - prevExonRank) > 1)
        {
            LOGGER.warn("transcript({}) invalid exon ranks(prev={} next={}) forwardStrand({}) position({})",
                    first.TransName, prevExonRank, nextExonRank, isForwardStrand, position);
            return null;
        }

        return new Transcript(geneAnnotation,
                first.TransName,
                isForwardStrand ? prevExonRank: nextExonRank, isForwardStrand ? prevExonEndPhase : nextExonEndPhase,
                isForwardStrand ? nextExonRank : prevExonEndPhase, isForwardStrand ? nextExonPhase : prevExonEndPhase,
                codingBases, totalCodingBases,
                exonMax, first.IsCanonical, transcriptStart, transcriptEnd,
                first.CodingStart, first.CodingEnd);
    }

    private List<HmfTranscriptRegion> findGeneRegions(long position, byte orientation, SortedSet<HmfTranscriptRegion> geneRegions)
    {
        List<HmfTranscriptRegion> matchedGenes = Lists.newArrayList();

        for(final HmfTranscriptRegion region : geneRegions)
        {
            if(position >= region.start() - PRE_GENE_PROMOTOR_DISTANCE && position <= region.end() + PRE_GENE_PROMOTOR_DISTANCE)
            {
                matchedGenes.add(region);
            }
        }

        return matchedGenes;
    }

    public static final String getSampleGeneAnnotationsFilename(final String path, final String sampleId)
    {
        String filename = path;

        if(!path.endsWith("/"))
                filename += File.separator;

        filename += sampleId + "_" + SV_GENE_TRANSCRIPTS_FILE_SUFFIX;

        return filename;
    }

    private static int VAR_ID_COL_INDEX = 0;
    private static int VAR_CHR_COL_INDEX = 1;
    private static int VAR_POS_COL_INDEX = 2;
    private static int VAR_ORIENT_COL_INDEX = 3;

    // gene data: isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
    private static int GENE_IS_START_COL_INDEX = 4;
    private static int GENE_NAME_COL_INDEX = 5;
    private static int GENE_STABLE_ID_COL_INDEX = 6;
    private static int GENE_STRAND_INDEX = 7;
    private static int GENE_SYNS_COL_INDEX = 8;
    private static int GENE_EIDS_COL_INDEX = 9;
    private static int GENE_KARYOTYPE_COL_INDEX = 10;

    // transcript data: transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, codingBase, totalCodingBases, exonMax, canonical, codingStart, codingEnd
    private static int TRANSCRIPT_ID_COL_INDEX = 11;
    private static int TRANSCRIPT_EUP_RANK_COL_INDEX = 12;
    private static int TRANSCRIPT_EUP_PHASE_COL_INDEX = 13;
    private static int TRANSCRIPT_EDN_RANK_COL_INDEX = 14;
    private static int TRANSCRIPT_EDN_PHASE_COL_INDEX = 15;
    private static int TRANSCRIPT_CDB_COL_INDEX = 16;
    private static int TRANSCRIPT_TCB_COL_INDEX = 17;
    private static int TRANSCRIPT_EMAX_COL_INDEX = 18;
    private static int TRANSCRIPT_CAN_COL_INDEX = 19;
    private static int TRANSCRIPT_TRANS_S_COL_INDEX = 20;
    private static int TRANSCRIPT_TRANS_E_COL_INDEX = 21;
    private static int TRANSCRIPT_CODE_S_COL_INDEX = 22;
    private static int TRANSCRIPT_CODE_E_COL_INDEX = 23;

    public boolean loadSampleGeneTranscripts(final String sampleId)
    {
        mSvIdGeneTranscriptsMap.clear();

        if(sampleId.isEmpty() || mDataPath.isEmpty())
            return false;

        final String filename = getSampleGeneAnnotationsFilename(mDataPath, sampleId);

        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty ensembl data file({})", filename);
                return false;
            }

            int currentVarId = -1;

            GeneAnnotation currentGene = null;
            List<GeneAnnotation> geneAnnotations = null;

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final int varId = Integer.parseInt(items[VAR_ID_COL_INDEX]);

                if(varId != currentVarId)
                {
                    if(currentVarId >= 0)
                    {
                        mSvIdGeneTranscriptsMap.put(currentVarId, geneAnnotations);
                    }

                    currentVarId = varId;
                    currentGene = null;

                    // start a new list for the new variant
                    geneAnnotations = Lists.newArrayList();
                }

                // isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                final String geneName = items[GENE_NAME_COL_INDEX];
                boolean geneIsStart = Boolean.parseBoolean(items[GENE_IS_START_COL_INDEX]);

                if(currentGene == null || !currentGene.geneName().equals(geneName) || currentGene.isStart() != geneIsStart)
                {
                    String[] synonymsStr = items[GENE_SYNS_COL_INDEX].split(";");
                    final List<String> synonyms = Lists.newArrayList(synonymsStr);

                    String[] entrezIdStr = items[GENE_EIDS_COL_INDEX].split(";");

                    final List<Integer> entrezIds = Lists.newArrayList();

                    for (int i = 0; i < entrezIdStr.length; ++i)
                    {
                        if(!entrezIdStr[i].isEmpty())
                            entrezIds.add(Integer.parseInt(entrezIdStr[i]));
                    }

                    currentGene = new GeneAnnotation(
                            varId,
                            geneIsStart,
                            geneName,
                            items[GENE_STABLE_ID_COL_INDEX],
                            Integer.parseInt(items[GENE_STRAND_INDEX]),
                            synonyms,
                            entrezIds,
                            items[GENE_KARYOTYPE_COL_INDEX]);

                    currentGene.setPositionalData(
                            items[VAR_CHR_COL_INDEX],
                            Long.parseLong(items[VAR_POS_COL_INDEX]),
                            Byte.parseByte(items[VAR_ORIENT_COL_INDEX]));

                    geneAnnotations.add(currentGene);
                }

                final String transcriptId = items[TRANSCRIPT_ID_COL_INDEX];


                int exonUpstreamRank = Integer.parseInt(items[TRANSCRIPT_EUP_RANK_COL_INDEX]);
                int exonUpstreamPhase = Integer.parseInt(items[TRANSCRIPT_EUP_PHASE_COL_INDEX]);
                int exonDownstreamRank = Integer.parseInt(items[TRANSCRIPT_EDN_RANK_COL_INDEX]);
                int exonDownstreamPhase = Integer.parseInt(items[TRANSCRIPT_EDN_PHASE_COL_INDEX]);

                // corrections for errors in Ensembl annotations

                if(exonDownstreamRank == -1 || exonUpstreamRank == -1 || abs(exonUpstreamRank - exonDownstreamRank) > 1)
                {
                    LOGGER.warn("skipping invalid transcript info: SV({}) trans({}) ranks(up={} down={})",
                            varId, transcriptId, exonUpstreamRank, exonDownstreamRank);
                }
                else
                {
                    // transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonMax, canonical, codingStart, codingEnd
                    Transcript transcript = new Transcript(
                            currentGene, transcriptId,
                            exonUpstreamRank, exonUpstreamPhase, exonDownstreamRank, exonDownstreamPhase,
                            Long.parseLong(items[TRANSCRIPT_CDB_COL_INDEX]),
                            Long.parseLong(items[TRANSCRIPT_TCB_COL_INDEX]),
                            Integer.parseInt(items[TRANSCRIPT_EMAX_COL_INDEX]),
                            Boolean.parseBoolean(items[TRANSCRIPT_CAN_COL_INDEX]),
                            Integer.parseInt(items[TRANSCRIPT_TRANS_S_COL_INDEX]),
                            Integer.parseInt(items[TRANSCRIPT_TRANS_E_COL_INDEX]),
                            items[TRANSCRIPT_CODE_S_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CODE_S_COL_INDEX]),
                            items[TRANSCRIPT_CODE_E_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CODE_E_COL_INDEX]));

                    currentGene.addTranscript(transcript);
                }

                line = fileReader.readLine();

                if(line == null)
                {
                    // add the last variant gene list
                    mSvIdGeneTranscriptsMap.put(varId, geneAnnotations);
                    break;
                }
            }

        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load sample gene annotations({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public void writeAnnotations(final String sampleId, final List<StructuralVariantAnnotation> annotations)
    {
        if(mDataPath.isEmpty() || sampleId.isEmpty())
            return;

        LOGGER.debug("writing {} annotations to file", annotations.size());

        String outputFilename = getSampleGeneAnnotationsFilename(mDataPath, sampleId);

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            // write header
            writer.write("SvId,Chromosome,Position,Orientation");
            writer.write(",IsStart,GeneName, GeneStableId, GeneStrand, Synonyms, EntrezIds, KaryotypeBand");
            writer.write(",TranscriptId,ExonUpstream,ExonUpstreamPhase,ExonDownstream,ExonDownstreamPhase,CodingBases,TotalCodingBases");
            writer.write(",ExonMax,Canonical,TranscriptStart,TranscriptEnd,CodingStart,CodingEnd,RegionType,CodingType");
            writer.newLine();

            for(final StructuralVariantAnnotation annotation : annotations)
            {
                if(annotation.annotations().isEmpty())
                {
                    // LOGGER.debug("SV({}) has no annotations", annotation.variant().primaryKey());
                    continue;
                }

                for(final GeneAnnotation geneAnnotation : annotation.annotations())
                {
                    String synonymnsStr = "";
                    for(final String syn : geneAnnotation.synonyms())
                    {
                        if(!synonymnsStr.isEmpty())
                            synonymnsStr += ";";

                        synonymnsStr += syn;
                    }

                    String entrezIdsStr = "";
                    for(final Integer eId : geneAnnotation.entrezIds())
                    {
                        if(!entrezIdsStr.isEmpty())
                            entrezIdsStr += ";";

                        entrezIdsStr += eId;
                    }

                    for(final Transcript transcript : geneAnnotation.transcripts())
                    {
                        final StructuralVariant var = annotation.variant();

                        boolean isStart = geneAnnotation.isStart();

                        writer.write(String.format("%d,%s,%d,%d",
                                var.primaryKey(), var.chromosome(isStart), var.position(isStart), var.orientation(isStart)));

                        // Gene info: isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                        writer.write(
                                String.format(",%s,%s,%s,%d,%s,%s,%s",
                                        geneAnnotation.isStart(),
                                        geneAnnotation.geneName(),
                                        geneAnnotation.stableId(),
                                        geneAnnotation.strand(),
                                        synonymnsStr,
                                        entrezIdsStr,
                                        geneAnnotation.karyotypeBand()));

                        // Transcript info: transcriptId,exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonStart, exonEnd, exonMax, canonical, codingStart, codingEnd
                        writer.write(
                                String.format(",%s,%d,%d,%d,%d,%d,%d",
                                        transcript.transcriptId(),
                                        transcript.exonUpstream(),
                                        transcript.exonUpstreamPhase(),
                                        transcript.exonDownstream(),
                                        transcript.exonDownstreamPhase(),
                                        transcript.codingBases(),
                                        transcript.totalCodingBases()));

                        writer.write(
                                String.format(",%d,%s,%d,%d,%d,%d,%s,%s",
                                        transcript.exonMax(),
                                        transcript.isCanonical(),
                                        transcript.transcriptStart(),
                                        transcript.transcriptEnd(),
                                        transcript.codingStart(),
                                        transcript.codingEnd(),
                                        transcript.regionType(),
                                        transcript.codingType()));

                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations: {}", e.toString());
        }
    }

    public final List<GeneAnnotation> updateAnnotationsByPosition(final EnrichedStructuralVariant var)
    {
        for (Map.Entry<Integer, List<GeneAnnotation>> entry : mSvIdGeneTranscriptsMap.entrySet())
        {
            // find transcript data by a position match, and then re-insert into the new map with the new ID
            final List<GeneAnnotation> geneList = entry.getValue();
            final GeneAnnotation gene = geneList.get(0);

            boolean matched = true;

            if (gene.isStart() && gene.chromosome().equals(var.chromosome(true))
            && gene.position() == var.position(true) && gene.orientation() == var.orientation(true))
            {
                matched = true;
            }
            else if (gene.isEnd() && gene.chromosome().equals(var.chromosome(false))
            && gene.position() == var.position(false) && gene.orientation() == var.orientation(false))
            {
                matched = true;
            }
            else
            {
                continue;
            }

            for (final GeneAnnotation annotation : geneList)
            {
                annotation.setVarId(var.primaryKey());
            }

            return geneList;
        }

        return null;
    }
}
