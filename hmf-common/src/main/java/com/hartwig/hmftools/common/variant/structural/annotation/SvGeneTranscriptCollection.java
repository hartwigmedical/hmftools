package com.hartwig.hmftools.common.variant.structural.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvGeneTranscriptCollection
{
    private String mDataPath;

    private Map<Integer, List<GeneAnnotation>> mSvIdGeneTranscriptsMap;

    public static String SV_GENE_TRANSCRIPTS_FILE_SUFFIX = "sv_ensembl_data.csv";

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public SvGeneTranscriptCollection()
    {
        mSvIdGeneTranscriptsMap = new HashMap();
    }

    public final Map<Integer, List<GeneAnnotation>> getSvIdGeneTranscriptsMap() { return mSvIdGeneTranscriptsMap; }

    public void setDataPath(final String dataPath)
    {
        mDataPath = dataPath;
    }

    private static final String getSampleGeneAnnotationsFilename(final String path, final String sampleId)
    {
        String filename = path;

        if(!path.endsWith("/"))
                filename += "/";

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
    private static int TRANSCRIPT_EUS_COL_INDEX = 12;
    private static int TRANSCRIPT_EUP_COL_INDEX = 13;
    private static int TRANSCRIPT_EDS_COL_INDEX = 14;
    private static int TRANSCRIPT_EDP_COL_INDEX = 15;
    private static int TRANSCRIPT_CDB_COL_INDEX = 16;
    private static int TRANSCRIPT_TCB_COL_INDEX = 17;
    private static int TRANSCRIPT_EMAX_COL_INDEX = 18;
    private static int TRANSCRIPT_CAN_COL_INDEX = 19;
    private static int TRANSCRIPT_CS_COL_INDEX = 20;
    private static int TRANSCRIPT_CE_COL_INDEX = 21;

    public boolean loadSampleGeneTranscripts(final String sampleId)
    {
        mSvIdGeneTranscriptsMap.clear();

        if(sampleId.isEmpty() || mDataPath.isEmpty())
            return false;

        final String filename = getSampleGeneAnnotationsFilename(mDataPath, sampleId);

        if (filename.isEmpty())
            return false;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty copy number CSV file({})", filename);
                return false;
            }

            int currentVarId = -1;

            GeneAnnotation currentGene = null;
            List<GeneAnnotation> geneAnnotations = null;

            int fileIndex = 0;

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
                        mSvIdGeneTranscriptsMap.put(varId, geneAnnotations);
                    }

                    currentVarId = varId;

                    // start a new list for the new variant
                    geneAnnotations = Lists.newArrayList();
                }

                // isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                final String geneName = items[GENE_NAME_COL_INDEX];

                if(currentGene == null || !currentGene.geneName().equals(geneName))
                {
                    if(currentGene != null)
                    {
                        // add to annotation and prepare a new one
                        geneAnnotations.add(currentGene);
                    }

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
                            Boolean.parseBoolean(items[GENE_IS_START_COL_INDEX]),
                            geneName,
                            items[GENE_STABLE_ID_COL_INDEX],
                            Integer.parseInt(items[GENE_STRAND_INDEX]),
                            synonyms,
                            entrezIds,
                            items[GENE_KARYOTYPE_COL_INDEX]);
                }

                final String transcriptId = items[TRANSCRIPT_ID_COL_INDEX];

                // transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonMax, canonical, codingStart, codingEnd
                Transcript transcript = new Transcript(
                        currentGene,
                        transcriptId,
                        Integer.parseInt(items[TRANSCRIPT_EUS_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EUP_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EDS_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EDP_COL_INDEX]),
                        Long.parseLong(items[TRANSCRIPT_CDB_COL_INDEX]),
                        Long.parseLong(items[TRANSCRIPT_TCB_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EMAX_COL_INDEX]),
                        Boolean.parseBoolean(items[TRANSCRIPT_CAN_COL_INDEX]),
                        items[TRANSCRIPT_CS_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CS_COL_INDEX]),
                        items[TRANSCRIPT_CE_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CE_COL_INDEX]));

                currentGene.addTranscript(transcript);

                line = fileReader.readLine();

                if(line == null)
                {
                    // add the last variant gene list
                    geneAnnotations.add(currentGene);
                    mSvIdGeneTranscriptsMap.put(varId, geneAnnotations);
                    break;
                }

                ++fileIndex;
            }

        }
        catch(IOException e)
        {
            LOGGER.error("failed to load sample gene annotations({}): {}", filename, e.toString());
        }

        return true;
    }

    public void setSvData(final List<StructuralVariantData> variants)
    {
        for(final StructuralVariantData var : variants)
        {
            List<GeneAnnotation> geneAnnotations = mSvIdGeneTranscriptsMap.get(Integer.parseInt(var.id()));

            if(geneAnnotations == null)
                continue;

            for(GeneAnnotation gene : geneAnnotations)
            {
                gene.setSvData(var);
            }
        }
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
            writer.write(",IsStart,GeneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand");
            writer.write(",TranscriptId,ExonUpstream,ExonUpstreamPhase,ExonDownstream,ExonDownstreamPhase,CodingBases,TotalCodingBases");
            writer.write(",ExonMax,Canonical,CodingStart,CodingEnd,RegionType,CodingType");
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
                                String.format(",%d,%s,%d,%d,%s,%s",
                                        transcript.exonMax(),
                                        transcript.isCanonical(),
                                        transcript.codingStart(),
                                        transcript.codingEnd(),
                                        transcript.getRegionType(),
                                        transcript.getCodingType()));

                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations");
        }
    }
}
