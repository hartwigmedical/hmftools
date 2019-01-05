package com.hartwig.hmftools.svannotation;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

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

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<Integer, List<GeneAnnotation>> mSvIdGeneTranscriptsMap;
    private Map<String, List<TranscriptExonData>> mGeneTransExonDataMap;
    private Map<String, List<EnsemblGeneData>> mEnsemblChrGeneDataMap;

    private BufferedWriter mBreakendWriter;

    // to get a wider range of candidate genes and filter by promotor distance later on
    public static int PRE_GENE_PROMOTOR_DISTANCE = 10000;

    public static String SV_GENE_TRANSCRIPTS_FILE_SUFFIX = "sv_ensembl_data.csv";

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mSvIdGeneTranscriptsMap = new HashMap();
        mGeneTransExonDataMap = new HashMap();
        mEnsemblChrGeneDataMap = new HashMap();
        mBreakendWriter = null;
    }

    public final Map<Integer, List<GeneAnnotation>> getSvIdGeneTranscriptsMap() { return mSvIdGeneTranscriptsMap; }
    public final Map<String, List<TranscriptExonData>> getGeneExonDataMap() { return mGeneTransExonDataMap; }

    public final EnsemblGeneData getGeneData(final String geneName)
    {
        for(Map.Entry<String, List<EnsemblGeneData>> entry : mEnsemblChrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData : entry.getValue())
            {
                if(geneData.GeneName.equals(geneName))
                    return geneData;
            }
        }

        return null;
    }

    public List<TranscriptExonData> getTransExonData(final String geneId)
    {
        return mGeneTransExonDataMap.get(geneId);
    }

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

            List<EnsemblGeneData> geneList = null;
            String currentChr = "";
            int geneCount = 0;

            while (line != null)
            {
                String[] items = line.split(",");

                final String geneId = items[GD_ID];
                final String chromosome = items[GD_CHR];

                EnsemblGeneData geneData = new EnsemblGeneData(
                        geneId, items[GD_NAME], chromosome, Byte.parseByte(items[GD_STRAND]),
                        Long.parseLong(items[GD_START]), Long.parseLong(items[GD_END]),
                        items[GD_ENTREZ], items[GD_BAND], items[GD_SYN]);

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;
                    geneList = mEnsemblChrGeneDataMap.get(chromosome);

                    if(geneList == null)
                    {
                        geneList = Lists.newArrayList();
                        mEnsemblChrGeneDataMap.put(chromosome, geneList);
                    }
                }

                geneData.setListIndex(geneList.size());
                geneList.add(geneData);
                ++geneCount;

                line = fileReader.readLine();
            }

            LOGGER.debug("loaded {} gene records", geneCount);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load Ensembl gene ({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    private static int SPECIFIC_VAR_ID = -1;
    // private static int SPECIFIC_VAR_ID = 5417525;

    public List<GeneAnnotation> findGeneAnnotationsBySv(int svId, boolean isStart, final String chromosome, long position, byte orientation)
    {
        List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        if(svId == SPECIFIC_VAR_ID)
        {
            LOGGER.debug("specific SV({})", svId);
        }

        final List<EnsemblGeneData> geneRegions = mEnsemblChrGeneDataMap.get(chromosome);

        final List<EnsemblGeneData> matchedGenes = findGeneRegions(position, geneRegions);

        // now look up relevant transcript and exon information
        for(final EnsemblGeneData geneData : matchedGenes)
        {
            final List<TranscriptExonData> transExonDataList = mGeneTransExonDataMap.get(geneData.GeneId);

            if (transExonDataList == null || transExonDataList.isEmpty())
                continue;

            GeneAnnotation currentGene = new GeneAnnotation(svId, isStart, geneData.GeneName, geneData.GeneId,
                    geneData.Strand, geneData.Synonyms, geneData.EntrezIds, geneData.KaryotypeBand);

            // collect up all the relevant exons for each unique transcript to analyse as a collection
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

                Transcript transcript = extractTranscriptExonData(transcriptExons, position, currentGene);

                if(transcript != null)
                {
                    currentGene.addTranscript(transcript);

                    // annotate with preceding gene info if the up distance isn't set
                    if(transcript.exonDistanceUp() == -1)
                    {
                        EnsemblGeneData precedingGene = findPrecedingGene(geneData, geneRegions);
                        if(precedingGene != null)
                        {
                            currentGene.setPrecedingGeneId(precedingGene.GeneId);
                            int preDistance = (int)abs(precedingGene.GeneEnd - position);
                            transcript.setExonDistances(preDistance, transcript.exonDistanceDown());
                        }
                    }
                }

                if(j == transExonDataList.size() - 1)
                    break;

                i = j - 1;
            }

            geneAnnotations.add(currentGene);
        }

        return geneAnnotations;
    }

    private List<EnsemblGeneData> findGeneRegions(long position, List<EnsemblGeneData> geneDataList)
    {
        List<EnsemblGeneData> matchedGenes = Lists.newArrayList();

        for(final EnsemblGeneData geneData : geneDataList)
        {
            long geneStartRange = geneData.Strand == 1 ? geneData.GeneStart - PRE_GENE_PROMOTOR_DISTANCE : geneData.GeneStart;
            long geneEndRange = geneData.Strand == 1 ? geneData.GeneEnd : geneData.GeneEnd + PRE_GENE_PROMOTOR_DISTANCE;

            if(position >= geneStartRange && position <= geneEndRange)
            {
                matchedGenes.add(geneData);
            }
        }

        return matchedGenes;
    }

    private EnsemblGeneData findPrecedingGene(final EnsemblGeneData geneData, List<EnsemblGeneData> geneDataList)
    {
        if(geneData.Strand == 1)
        {
            for(int i = geneData.getListIndex() - 1; i >= 0; --i)
            {
                final EnsemblGeneData gene = geneDataList.get(i);
                if(gene.Strand == geneData.Strand)
                    return gene;
            }
        }
        else
        {
            for(int i = geneData.getListIndex() + 1; i < geneDataList.size(); ++i)
            {
                final EnsemblGeneData gene = geneDataList.get(i);
                if(gene.Strand == geneData.Strand)
                    return gene;
            }
        }

        return null;
    }

    private Transcript extractTranscriptExonData(final List<TranscriptExonData> transcriptExons, long position, final GeneAnnotation geneAnnotation)
    {
        int exonMax = transcriptExons.size();

        final TranscriptExonData first = transcriptExons.get(0);

        boolean isForwardStrand = geneAnnotation.strand() == 1;

        int upExonRank = -1;
        int upExonPhase = -1;
        int downExonRank = -1;
        int downExonPhase = -1;
        long nextUpDistance = -1;
        long nextDownDistance = -1;

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

        for (int index = 0; index < transcriptExons.size(); ++index)
        {
            final TranscriptExonData exonData = transcriptExons.get(index);
            long exonStart = exonData.ExonStart;
            long exonEnd = exonData.ExonEnd;

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

        if(!isForwardStrand)
        {
            codingBases = totalCodingBases - codingBases;
        }

        Transcript transcript = new Transcript(geneAnnotation,
                first.TransName,
                upExonRank, upExonPhase, downExonRank, downExonPhase,
                codingBases, totalCodingBases,
                exonMax, first.IsCanonical, first.TransStart, first.TransEnd,
                first.CodingStart, first.CodingEnd);

        transcript.setBioType(first.BioType);
        transcript.setExonDistances((int)nextUpDistance, (int)nextDownDistance);

        return transcript;
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

                // write header
                mBreakendWriter.write("SampleId,SvId,IsStart,Chromosome,Position,Orientation,Type");
                mBreakendWriter.write(",GeneName,GeneStableId,GeneStrand,TranscriptId,IsCanonical,BioType,TransStart,TransEnd");
                mBreakendWriter.write(",ExonRankUp,ExonPhaseUp,ExonRankDown,ExonPhaseDown,CodingBases,TotalCodingBases");
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
                                        geneAnnotation.geneName(), geneAnnotation.stableId(), geneAnnotation.strand(),
                                        transcript.transcriptId(), transcript.isCanonical(), transcript.bioType(),
                                        transcript.transcriptStart(), transcript.transcriptEnd()));

                        // Transcript info: exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonStart, exonEnd, exonMax, canonical, codingStart, codingEnd
                        writer.write(
                                String.format(",%d,%d,%d,%d,%d,%d",
                                        transcript.exonUpstream(), transcript.exonUpstreamPhase(),
                                        transcript.exonDownstream(), transcript.exonDownstreamPhase(),
                                        transcript.codingBases(), transcript.totalCodingBases()));

                        writer.write(
                                String.format(",%d,%d,%d,%s,%s,%d,%d",
                                        transcript.exonMax(), transcript.codingStart(), transcript.codingEnd(),
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

}
