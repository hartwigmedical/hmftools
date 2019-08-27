package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Svbreakend;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ensembl.database.homo_sapiens_core.tables.Gene;

public class DisruptionFinder
{
    private final SvGeneTranscriptCollection mGeneTransCollection;
    private Set<String> mDisruptionGeneIDPanel;

    private final List<Transcript> mDisruptions;

    private boolean mNewDisruptionLogic;

    private BufferedWriter mWriter;
    private final String mOutputDir;

    public static final String DRUP_TSG_GENES_FILE = "drup_tsg_file";

    private static final Logger LOGGER = LogManager.getLogger(DisruptionFinder.class);

    public DisruptionFinder(final CommandLine cmd, final SvGeneTranscriptCollection geneTransCache, final String outputDir)
    {
        mGeneTransCollection = geneTransCache;
        mDisruptionGeneIDPanel = null;

        mDisruptions = Lists.newArrayList();

        mNewDisruptionLogic = false;

        mOutputDir = outputDir;
        mWriter = null;

        initialise(cmd);
    }

    public final List<Transcript> getDisruptions() { return mDisruptions; }

    public void setNewDisruptionLogic(boolean toggle) { mNewDisruptionLogic = toggle; }

    private void initialise(final CommandLine cmd)
    {
        mDisruptionGeneIDPanel = tsgDriverGeneIDs();

        // TEMP: load DRUP TSGs from file
        if(cmd != null && cmd.hasOption(DRUP_TSG_GENES_FILE))
        {
            loadDrupTSGs(cmd.getOptionValue(DRUP_TSG_GENES_FILE));
        }
    }

    public boolean matchesDisruptionGene(final GeneAnnotation gene)
    {
        return mDisruptionGeneIDPanel.stream().anyMatch(geneID -> gene.synonyms().contains(geneID));
    }

    public void addDisruptionGene(final String geneId)
    {
        if(!mDisruptionGeneIDPanel.contains(geneId))
            mDisruptionGeneIDPanel.add(geneId);
    }

    public void findReportableDisruptions(final List<SvVarData> svList)
    {
        mDisruptions.clear();

        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                if (be == SE_END && var.isSglBreakend())
                    continue;

                final List<GeneAnnotation> tsgGenesList = var.getGenesList(isStart(be)).stream()
                        .filter(x -> matchesDisruptionGene(x)).collect(Collectors.toList());

                for(GeneAnnotation gene : tsgGenesList)
                {
                    List<Transcript> reportableDisruptions = gene.transcripts().stream()
                            .filter(Transcript::isCanonical)
                            .filter(Transcript::isDisruptive)
                            .collect(Collectors.toList());

                    for(Transcript transcript : reportableDisruptions)
                    {
                        LOGGER.debug("var({}) breakend({}) gene({}) transcript({}) is disrupted",
                                var.id(), var.getBreakend(be), gene.GeneName, transcript.StableId);

                        transcript.setReportableDisruption(true);
                        mDisruptions.add(transcript);
                    }
                }
            }
        }
    }

    // to be deprecated
    public void markNonDisruptiveTranscripts(final SvVarData var)
    {
        if(!var.isSimpleType())
            return;

        final List<GeneAnnotation> genesStart = var.getGenesList(true);
        final List<GeneAnnotation> genesEnd = var.getGenesList(false);

        if(genesStart.isEmpty() || genesEnd.isEmpty())
            return;

        for(final GeneAnnotation geneStart : genesStart)
        {
            final GeneAnnotation geneEnd = genesEnd.stream()
                    .filter(x -> x.StableId.equals(geneStart.StableId)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            for (final Transcript transStart : geneStart.transcripts())
            {
                final Transcript transEnd = geneEnd.transcripts().stream()
                        .filter(x -> x.StableId.equals(transStart.StableId)).findFirst().orElse(null);

                if(transEnd == null)
                    continue;

                if(transStart.ExonUpstream == transEnd.ExonUpstream)
                {
                    transStart.setIsDisruptive(false);
                    transEnd.setIsDisruptive(false);
                }
            }
        }
    }

    public void markTranscriptsDisruptive(final SvVarData var)
    {
        if(!mNewDisruptionLogic)
        {
            markNonDisruptiveTranscripts(var);
            return;
        }

        if(var.isSglBreakend())
            return;

        final List<GeneAnnotation> genesStart = var.getGenesList(true);
        final List<GeneAnnotation> genesEnd = var.getGenesList(false);

        if(genesStart.isEmpty() && genesEnd.isEmpty())
            return;

        /* test the breakend to see if:
            - it isn't chained - revert to simple disrupted definitions
            - it is chained
                - the chain returns to the same intro
                - the chain traverses a splice acceptor in any direction
        */

        final SvCluster cluster = var.getCluster();

        boolean isSimpleSV = var.isSimpleType();

        if(isSimpleSV)
        {
            markNonDisruptiveGeneTranscripts(genesStart, genesEnd);
        }

        final List<SvChain> chains = cluster.findChains(var);

        // first test each breakend in turn to see if it forms a TI wholy within an intron and where the other breakends
        // of the TI are non-genic

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final SvBreakend breakend = var.getBreakend(se);

            final List<GeneAnnotation> svGenes = se == SE_START ? genesStart : genesEnd;

            List<Transcript> transList = getDisruptedTranscripts(svGenes);

            if(transList.isEmpty())
                continue;

            boolean otherBreakendNonGenic = se == SE_START ? genesEnd.isEmpty() : genesStart.isEmpty();

            final List<SvLinkedPair> links = var.getLinkedPairs(breakend.usesStart());

            for(final SvLinkedPair pair : links)
            {
                final SvBreakend otherSvBreakend = pair.getOtherBreakend(breakend);

                List<GeneAnnotation> otherGenes = otherSvBreakend.getSV().getGenesList(otherSvBreakend.usesStart());
                List<Transcript> otherTransList = getDisruptedTranscripts(otherGenes);

                boolean otherSvOtherBreakendNonGenic = otherSvBreakend.getSV().getGenesList(!otherSvBreakend.usesStart()).isEmpty();

                if(otherBreakendNonGenic && otherSvOtherBreakendNonGenic)
                {
                    markNonDisruptiveTranscripts(transList, otherTransList);
                    removeNonDisruptedTranscripts(transList);
                }
            }

            // next test for a breakend which returns to the same intro via a chain with the same orientation
            // without passing through any splice acceptors
            if(!transList.isEmpty())
            {
                for(SvChain chain : chains)
                {
                    checkChainedTranscripts(breakend, chain, transList);
                }
            }
        }
    }

    private static final List<Transcript> getDisruptedTranscripts(final List<GeneAnnotation> genes)
    {
        List<Transcript> transList = Lists.newArrayList();

        for (final GeneAnnotation gene : genes)
        {
            transList.addAll(gene.transcripts().stream().filter(x -> x.isDisruptive()).collect(Collectors.toList()));
        }

        return transList;
    }

    private static void removeNonDisruptedTranscripts(final List<Transcript> transList)
    {
        int index = 0;

        while(index < transList.size())
        {
            if(transList.get(index).isDisruptive())
                ++index;
            else
                transList.remove(index);
        }
    }

    private void checkChainedTranscripts(final SvBreakend breakend, final SvChain chain, final List<Transcript> transList)
    {
        // an SV whose breakends are not both within the same intron disrupts those transcripts unless
        // a) one breakend isn't genic and the other forms a TI whole within an intron OR
        // b) both breakends are in chained sections which come back to the same intro with correct orientation and
        // without traversing a splice acceptor OR

        final List<SvLinkedPair> links = chain.getLinkedPairs();

        for(int i = 0; i < links.size(); ++i)
        {
            int index = 0;
            boolean traverseUp = false;

            if(i == 0 && chain.getOpenBreakend(true) == breakend)
            {
                index = -1;
                traverseUp = true;
            }
            else if(i == links.size() - 1 && chain.getOpenBreakend(false) == breakend)
            {
                index = links.size();
                traverseUp = false;
            }
            else if(links.get(i).hasBreakend(breakend))
            {
                index = i;
                traverseUp = links.get(i).secondBreakend() == breakend;
            }
            else
            {
                continue;
            }

            while(true)
            {
                index += traverseUp ? 1 : -1;

                if(index < 0 || index >= links.size())
                    break;

                final SvLinkedPair nextPair = links.get(index);

                // does this next link traverse another splice acceptor?
                boolean traversesGene = pairTraversesGene(nextPair, 0, false);

                if(traversesGene)
                    return;

                // no need to chain this breakend any further as soon as any of its chains have crossed another splice acceptor

                // does it return to the same intron and correct orientation for any transcripts
                final SvBreakend nextBreakend = traverseUp ?
                        nextPair.secondBreakend().getOtherBreakend() : nextPair.firstBreakend().getOtherBreakend();

                if(nextBreakend.orientation() == breakend.orientation() || !nextBreakend.chromosome().equals(breakend.chromosome()))
                    continue;

                List<GeneAnnotation> otherGenes = nextBreakend.getSV().getGenesList(nextBreakend.usesStart());

                if(otherGenes.isEmpty())
                    continue;

                List<Transcript> otherTransList = getDisruptedTranscripts(otherGenes);

                if(!otherTransList.isEmpty())
                {
                    markNonDisruptiveTranscripts(transList, otherTransList);
                    removeNonDisruptedTranscripts(transList);

                    if(transList.isEmpty())
                        return;
                }
            }
        }
    }

    private static void markNonDisruptiveGeneTranscripts(final List<GeneAnnotation> genesStart, final List<GeneAnnotation> genesEnd)
    {
        for(final GeneAnnotation geneStart : genesStart)
        {
            final GeneAnnotation geneEnd = genesEnd.stream()
                    .filter(x -> x.StableId.equals(geneStart.StableId)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            markNonDisruptiveTranscripts(geneStart.transcripts(), geneEnd.transcripts());
        }
    }

    private static void markNonDisruptiveTranscripts(final List<Transcript> transList1, final List<Transcript> transList2)
    {
        for (final Transcript trans1 : transList1)
        {
            final Transcript transEnd = transList2.stream()
                    .filter(x -> x.StableId.equals(trans1.StableId)).findFirst().orElse(null);

            if(transEnd == null)
                continue;

            if(trans1.ExonUpstream == transEnd.ExonUpstream)
            {
                trans1.setIsDisruptive(false);
                transEnd.setIsDisruptive(false);
            }
        }
    }

    public boolean pairTraversesGene(SvLinkedPair pair, int fusionDirection, boolean isPrecodingUpstream)
    {
        // for this pair to not affect the fusion, the section it traverses cannot cross any gene's splice acceptor
        // with the same strand direction unless that is part of a fully traversed non-coding 5' exon

        long lowerPos = pair.getBreakend(true).position();
        long upperPos = pair.getBreakend(false).position();

        List<EnsemblGeneData> geneDataList = mGeneTransCollection.getChrGeneDataMap().get(pair.chromosome());

        if(geneDataList == null)
            return false;

        for(EnsemblGeneData geneData : geneDataList)
        {
            if(lowerPos > geneData.GeneEnd)
                continue;

            if(upperPos < geneData.GeneStart)
                break;

            if(fusionDirection == 0 || geneData.Strand == fusionDirection)
            {
                // check whether a splice acceptor is encountered within this window
                List<TranscriptData> transDataList = mGeneTransCollection.getTranscripts(geneData.GeneId);

                if(transDataList == null)
                    continue;

                for(final TranscriptData transData : transDataList)
                {
                    for (final ExonData exonData : transData.exons())
                    {
                        if (exonData.ExonRank == 1)
                            continue;

                        if ((geneData.Strand == 1 && lowerPos <= exonData.ExonStart && upperPos >= exonData.ExonStart)
                                || (geneData.Strand == -1 && lowerPos <= exonData.ExonEnd && upperPos >= exonData.ExonEnd))
                        {
                            // allow an exon to be fully traversed if the upstream transcript is pre-coding
                            if (isPrecodingUpstream && lowerPos <= exonData.ExonStart && upperPos >= exonData.ExonEnd)
                            {
                                if (geneData.Strand == 1 && (transData.CodingStart == null || upperPos < transData.CodingStart))
                                    continue;
                                else if (geneData.Strand == -1 && (transData.CodingEnd == null || lowerPos > transData.CodingEnd))
                                    continue;
                            }

                            LOGGER.trace("pair({}) direction({}) traverses splice acceptor({} {}) exon(rank{} pos={})",
                                    pair.toString(), fusionDirection, geneData.GeneName, transData.TransName,
                                    exonData.ExonRank, exonData.ExonStart, exonData.ExonEnd);

                            return true;
                        }
                    }
                }
            }
        }

        return false;
    }

    private static boolean areDisruptivePair(final Transcript trans1, final Transcript trans2)
    {
        if(trans1.parent().id() != trans2.parent().id())
            return true;

        if(!trans1.StableId.equals(trans2.StableId))
            return true;

        // only DELs, DUPs and INS
        if(trans1.parent().orientation() == trans2.parent().orientation())
            return true;

        if(trans1.ExonUpstream != trans2.ExonUpstream)
            return true;

        return false;
    }

    private static Set<String> tsgDriverGeneIDs()
    {
        Set<String> tsgDriverGeneIDs = Sets.newHashSet();
        Map<String, HmfTranscriptRegion> allGenes = HmfGenePanelSupplier.allGenesMap37();

        for (String gene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet())
        {
            tsgDriverGeneIDs.add(allGenes.get(gene).geneID());
        }

        return tsgDriverGeneIDs;
    }

    private void loadDrupTSGs(final String filename)
    {
        if (filename.isEmpty() || !Files.exists(Paths.get(filename)))
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty DRUP TSG file({})", filename);
                return;
            }

            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                // parse CSV data
                final String[] items = line.split(",");

                if(items.length < 2)
                {
                    LOGGER.error("invalid DRUP TSG record: {}", line);
                    return;
                }

                final String geneName = items[0];

                final EnsemblGeneData geneData = mGeneTransCollection.getGeneDataByName(geneName);
                if(geneData != null)
                {
                    addDisruptionGene(geneData.GeneId);
                }
                else
                {
                    LOGGER.error("gene data not found for gene({})", geneName);
                }

                line = fileReader.readLine();
            }
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to load DRUP TSG file({}): {}", filename, e.toString());
            return;
        }
    }

    public void writeSampleData(final String sampleId)
    {
        // write sample file for patient reporter
        List<ReportableDisruption> reportedDisruptions = Lists.newArrayList();

        for(final Transcript transcript : mDisruptions)
        {
            final GeneAnnotation gene = transcript.parent();

            reportedDisruptions.add(ImmutableReportableDisruption.builder()
                    .svId(gene.id())
                    .chromosome(gene.chromosome())
                    .orientation(gene.orientation())
                    .strand(gene.Strand)
                    .chrBand(gene.karyotypeBand())
                    .gene(transcript.geneName())
                    .type(gene.type().toString())
                    .ploidy(gene.ploidy())
                    .exonUp(transcript.ExonUpstream)
                    .exonDown(transcript.ExonDownstream)
                    .build());
        }

        try
        {
            final String disruptionsFile = ReportableDisruptionFile.generateFilename(mOutputDir, sampleId);
            ReportableDisruptionFile.write(disruptionsFile, reportedDisruptions);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write sample disruptions file: {}", e.toString());
        }
    }

    public void initialiseOutputFile(final String fileName)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = mOutputDir + fileName;

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,Reportable,SvId,IsStart,Chromosome,Position,Orientation");
                mWriter.write(",GeneId,GeneName,Strand,TransId,ExonUp,ExonDown,CodingType,RegionType");
                mWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing disruptions: {}", e.toString());
        }
    }

    public void writeMultiSampleData(final String sampleId)
    {
        if(mWriter == null)
            return;

        try
        {
            for(final Transcript transcript : mDisruptions)
            {
                final GeneAnnotation gene = transcript.parent();

                mWriter.write(String.format("%s,%s,%d,%s,%s,%d,%d",
                        sampleId, transcript.reportableDisruption(), gene.id(), gene.isStart(),
                        gene.chromosome(), gene.position(), gene.orientation()));

                mWriter.write(String.format(",%s,%s,%d,%s,%d,%d,%s,%s",
                        gene.StableId, gene.GeneName, gene.Strand, transcript.StableId,
                        transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(), transcript.regionType()));
            }

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }


}
