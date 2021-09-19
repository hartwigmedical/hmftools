package com.hartwig.hmftools.svtools.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.svtools.germline.GermlineVcfConfig.GENE_PANEL_FILE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.StructuralVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneImpact
{
    // private final GermlineVcfConfig mConfig;
    // private final EnsemblDataCache mGeneDataCache;

    private BufferedWriter mWriter;

    public static final String DISRUPTION_TYPE_NONE = "NONE";
    public static final String DISRUPTION_TYPE_OVERLAP = "OVERLAP";
    public static final String DISRUPTION_TYPE_ONE_BREAK = "ONE_BREAK";
    public static final String DISRUPTION_TYPE_TWO_BREAK_COMPLEX = "TWO_BREAK_COMPLEX";
    public static final String DISRUPTION_TYPE_TWO_BREAK_SIMPLE = "TWO_BREAK_SIMPLE";
    public static final String DISRUPTION_TYPE_RECIP_INV = "RECIP_INV";
    public static final String DISRUPTION_TYPE_SHARD = "SHARD";

    private static final Logger LOGGER = LogManager.getLogger(GeneImpact.class);

    /*

    public GeneImpact(final GermlineVcfConfig config, final CommandLine cmd)
    {
        mConfig = config;

        if(cmd.hasOption(ENSEMBL_DATA_DIR))
        {
            final List<String> genePanelIds = cmd.hasOption(GENE_PANEL_FILE) ?
                    loadGenePanel(cmd.getOptionValue(GENE_PANEL_FILE)) : Lists.newArrayList();

            mGeneDataCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), RefGenomeVersion.RG_37);
            mGeneDataCache.setRestrictedGeneIdList(genePanelIds);
            mGeneDataCache.load(true);

            mGeneDataCache.setRequiredData(config.CheckDisruptions, false, false, true);

            if(config.CheckDisruptions)
            {
                mGeneDataCache.loadTranscriptData(genePanelIds);
            }
        }
        else
        {
            mGeneDataCache = null;
        }

        mWriter = null;
    }

    public void close() { closeBufferedWriter(mWriter); }

    private final List<String> loadGenePanel(final String geneFile)
    {
        final List<String> genePanelIds = Lists.newArrayList();

        if (!Files.exists(Paths.get(geneFile)))
            return genePanelIds;

        try
        {
            List<String> fileContents = Files.readAllLines(new File(geneFile).toPath());

            for(final String geneData : fileContents)
            {
                if(geneData.contains("GeneId"))
                    continue;

                String[] items = geneData.split(",");
                if(items.length != 2)
                {
                    LOGGER.error("invalid geneData({}) - expected 'GeneId,GeneName'");
                    return genePanelIds;
                }

                genePanelIds.add(items[0]);
            }

            LOGGER.info("loaded genePanelFile({}) with {} genes", geneFile, fileContents.size());
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load gene panel file({}): {}", geneFile, e.toString());
        }

        return genePanelIds;
    }

    public static final int OVERLAPPED_GENES = 2;

    public void populateGeneAnnotations(
            final StructuralVariant sv, int svIndex, List<List<GeneAnnotation>> breakendPairGenes, List<GeneAnnotation> overlapGenes)
    {
        if(mGeneDataCache == null)
            return;

        // find any gene overlapped by a DEL or with a breakend in it
        if(sv.type() == DEL)
        {
            overlapGenes.addAll(mGeneDataCache.findGeneAnnotationsByOverlap(
                    svIndex, sv.chromosome(true), sv.position(true).intValue(), sv.position(false).intValue()));

            if(!overlapGenes.isEmpty())
            {
                overlapGenes.stream().forEach(x -> x.setPositionalData(
                        sv.chromosome(true), sv.position(true).intValue(), sv.orientation(true)));
            }
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(sv.type() == SGL && se == SE_END)
                continue;

            boolean isStart = isStart(se);

            List<GeneAnnotation> svGenes = breakendPairGenes.get(se);

            svGenes.addAll(mGeneDataCache.findGeneAnnotationsBySv(
                    svIndex, isStart, sv.chromosome(isStart), sv.position(isStart).intValue(), sv.orientation(isStart), 0)
                    .stream().filter(x -> x.canonical() != null).collect(Collectors.toList()));

            if(!svGenes.isEmpty())
            {
                svGenes.stream().forEach(x -> x.setPositionalData(
                        sv.chromosome(isStart), sv.position(isStart).intValue(), sv.orientation(isStart)));
            }
        }
    }

    public void findDisruptiveVariants(final String sampleId, final List<GermlineSV> germlineSVs)
    {
        // SVs are not disruptive if they:
        // - are outside a gene unless an overlapping DEL
        // - form a TI within a single intron and the other breakend is outside a gene
        // - are a DUP or DEL and contained within an intron

        for(GermlineSV var : germlineSVs)
        {
            final StructuralVariant sv = var.sv();

            final Map<GeneAnnotation,String> svDisruptions = var.getDisruptions();

            var.getOverlappedGenes().forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_OVERLAP));

            final List<GeneAnnotation> genesStart = var.getBreakendGenes(SE_START);
            final List<GeneAnnotation> genesEnd = var.getBreakendGenes(SE_END);

            if(genesStart.isEmpty() && genesEnd.isEmpty())
                continue;

            if(genesStart.isEmpty() || genesEnd.isEmpty())
            {
                final List<GeneAnnotation> genesList = !genesStart.isEmpty() ? genesStart : genesEnd;

                // one breakend disrupts a gene
                genesList.forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_ONE_BREAK));

                final Transcript trans = genesList.get(0).canonical();

                LOGGER.debug("var({}:{} {}:{} -> {}:{}) gene({}) disrupted by {} breakend exons({} & {})",
                        sv.id(), sv.type(), sv.chromosome(true), sv.position(true),
                        sv.end() != null ? sv.chromosome(false) : "0", sv.end() != null ? sv.position(false) : -1,
                        trans.geneName(), !genesList.isEmpty() ? "start" : "end", trans.ExonUpstream, trans.ExonDownstream);
            }
            else if(var.Type == BND || var.Type == INV)
            {
                // both ends are in genes and initially are considered disruptive
                genesStart.forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_TWO_BREAK_COMPLEX));
                genesEnd.forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_TWO_BREAK_COMPLEX));

                final Transcript transStart = genesStart.get(0).canonical();
                final Transcript transEnd = genesEnd.get(0).canonical();

                LOGGER.debug("var({}:{} {}:{} -> {}:{}) genes({} & {}) disrupted",
                        sv.id(), sv.type(), sv.chromosome(true), sv.position(true),
                        sv.end() != null ? sv.chromosome(false) : "0", sv.end() != null ? sv.position(false) : -1,
                        transStart != null ? transStart.geneName() : "none", transEnd != null ? transEnd.geneName() : "none");
            }
            else if(isDisruptiveSV(sv, genesStart, genesEnd))
            {
                genesStart.forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_TWO_BREAK_SIMPLE));
                genesEnd.forEach(x -> svDisruptions.put(x, DISRUPTION_TYPE_TWO_BREAK_SIMPLE));
                continue;
            }
        }

        // check for SVs linked within the same intron
        removeIntronicLinks(germlineSVs);

        // write output to file
        germlineSVs.stream().filter(x -> !x.getDisruptions().isEmpty()).forEach(x -> writeDisruptions(sampleId, x));
    }

    private boolean isDisruptiveSV(final StructuralVariant var, final List<GeneAnnotation> genesStart, final List<GeneAnnotation> genesEnd)
    {
        for(final GeneAnnotation geneStart : genesStart)
        {
            final GeneAnnotation geneEnd = genesEnd.stream()
                    .filter(x -> x.StableId.equals(geneStart.StableId)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            Transcript transStart = geneStart.canonical();
            Transcript transEnd = geneEnd.canonical();

            if(transStart == null || transEnd == null)
                continue;

            if(var.type() == DUP)
            {
                if(transStart.ExonUpstream == 1 && transEnd.ExonDownstream <= 2 && !transStart.isExonic())
                {
                    return false;
                }
            }

            if(var.type() == DUP || var.type() == DEL)
            {
                if (transStart.ExonUpstream == transEnd.ExonUpstream && !transStart.isExonic() && !transEnd.isExonic())
                {
                    return false;
                }
            }
            else if(var.type() == INS && !transStart.isExonic())
            {
                return false;
            }

            LOGGER.debug("var({}:{} {}:{} -> {}:{}) gene({}) disrupted in exons({} & {})",
                    var.id(), var.type(), var.chromosome(true), var.position(true),
                    var.end() != null ? var.chromosome(false) : "0", var.end() != null ? var.position(false) : -1,
                    geneStart.GeneName, transStart.ExonUpstream, transEnd.ExonUpstream);

            return true;
        }

        return false;
    }

    private void removeIntronicLinks(final List<GermlineSV> germlineSVs)
    {
        final List<GermlineSV> recipInvSVs = Lists.newArrayList();
        final List<GermlineSV> intronicShardSVs = Lists.newArrayList();

        for(int i = 0; i < germlineSVs.size() - 1; ++i)
        {
            final GermlineSV var1 = germlineSVs.get(i);

            if(intronicShardSVs.contains(var1) || recipInvSVs.contains(var1))
                continue;

            final List<GeneAnnotation> genesStart1 = var1.getBreakendGenes(SE_START);
            final List<GeneAnnotation> genesEnd1 = var1.getBreakendGenes(SE_END);

            boolean dualBreakends1 = !genesStart1.isEmpty() && !genesEnd1.isEmpty();

            if(dualBreakends1 && var1.Type != INV)
                continue;

            final List<GeneAnnotation> genesList1 = Lists.newLinkedList(genesStart1);
            genesList1.addAll(genesEnd1);

            for(int j = i + 1; j < germlineSVs.size(); ++j)
            {
                final GermlineSV var2 = germlineSVs.get(j);

                if(intronicShardSVs.contains(var2) || recipInvSVs.contains(var2))
                    continue;

                final List<GeneAnnotation> genesStart2 = var2.getBreakendGenes(SE_START);
                final List<GeneAnnotation> genesEnd2 = var2.getBreakendGenes(SE_END);

                boolean dualBreakends2 = !genesStart2.isEmpty() && !genesEnd2.isEmpty();

                if(dualBreakends2 && var2.Type != INV)
                    continue;

                final List<GeneAnnotation> genesList2 = Lists.newLinkedList(genesStart2);
                genesList2.addAll(genesEnd2);

                // first check if the 2 breakend transcripts are facing within the same intron
                if(!dualBreakends1 && !dualBreakends2)
                {
                    for (final GeneAnnotation gene1 : genesList1)
                    {
                        final Transcript trans1 = gene1.canonical();
                        final GeneAnnotation gene2 =
                                genesList2.stream().filter(x -> x.StableId.equals(gene1.StableId)).findFirst().orElse(null);

                        if (gene2 == null)
                            continue;

                        final Transcript trans2 = gene2.canonical();

                        if (!trans1.isIntronic() || !trans2.isIntronic())
                            continue;

                        if (trans1.ExonDownstream != trans2.ExonDownstream)
                            continue;

                        if (gene1.orientation() == gene2.orientation())
                            continue;

                        if ((gene1.orientation() == -1 && gene1.position() < gene2.position())
                        || (gene2.orientation() == -1 && gene2.position() < gene1.position()))
                        {
                            int tiLength = abs(gene2.position() - gene1.position());

                            LOGGER.debug("SVs({} & {}) have facing intronic breakends in gene({}) exons({} -> {}) tiLength({})",
                                    var1.Id, var2.Id, gene1.GeneName, trans1.ExonUpstream, trans1.ExonDownstream, tiLength);

                            intronicShardSVs.add(var1);
                            intronicShardSVs.add(var2);

                            var1.getDisruptions().put(gene1, DISRUPTION_TYPE_SHARD);
                            var2.getDisruptions().put(gene2, DISRUPTION_TYPE_SHARD);
                        }
                    }
                }
                else if(dualBreakends1 && dualBreakends2)
                {
                    for (final GeneAnnotation gene1 : genesStart1)
                    {
                        final GeneAnnotation geneEnd1 = genesEnd1.stream()
                                .filter(x -> x.StableId.equals(gene1.StableId)).findFirst().orElse(null);

                        if (geneEnd1 == null)
                            continue;

                        final GeneAnnotation gene2 =
                                genesStart2.stream().filter(x -> x.StableId.equals(gene1.StableId)).findFirst().orElse(null);

                        final GeneAnnotation geneEnd2 =
                                genesEnd2.stream().filter(x -> x.StableId.equals(gene1.StableId)).findFirst().orElse(null);

                        if (gene2 == null || geneEnd2 == null)
                            continue;

                        final Transcript transEnd1 = geneEnd1.canonical();
                        final Transcript transStart1 = gene1.canonical();
                        final Transcript transStart2 = gene2.canonical();
                        final Transcript transEnd2 = geneEnd2.canonical();

                        if (!transStart1.isIntronic() || !transEnd1.isIntronic() || !transStart2.isIntronic() || !transEnd2.isIntronic())
                            continue;

                        // check orientations are opposite
                        if (transStart1.gene().orientation() == transStart2.gene().orientation())
                            continue;

                        if (transStart1.geneName().equals(transEnd1.geneName()) && transStart1.ExonDownstream == transEnd1.ExonDownstream
                                && transStart1.geneName().equals(transStart2.geneName())
                                && transStart1.ExonDownstream == transStart2.ExonDownstream
                                && transStart1.geneName().equals(transEnd2.geneName())
                                && transStart1.ExonDownstream == transEnd2.ExonDownstream)
                        {
                            int dbLength = abs(transStart1.gene().position() - transStart2.gene().position());

                            LOGGER.debug("SVs({} & {}) form reciprocal INV in gene({}) exons({} -> {}) dbLength({})",
                                    var1.Id, var2.Id, transStart1.geneName(),
                                    transStart1.ExonUpstream, transStart1.ExonDownstream, dbLength);

                            var1.getDisruptions().put(gene1, DISRUPTION_TYPE_RECIP_INV);
                            var1.getDisruptions().put(geneEnd1, DISRUPTION_TYPE_RECIP_INV);
                            var2.getDisruptions().put(gene2, DISRUPTION_TYPE_RECIP_INV);
                            var2.getDisruptions().put(geneEnd2, DISRUPTION_TYPE_RECIP_INV);

                            recipInvSVs.add(var1);
                            recipInvSVs.add(var2);
                        }
                    }
                }
            }
        }
    }

    private void writeDisruptions(final String sampleId, final GermlineSV germlineSV)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = mConfig.OutputDir + "LNX_GERMLINE_DISRUPTIONS.csv";

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,SvId,IsStart,Type,Chromosome,Position,Orientation");
                mWriter.write(",GeneId,GeneName,TransId,ExonUp,ExonDown,CodingType,RegionType,DisruptionType");
                mWriter.newLine();
            }

            final Map<GeneAnnotation,String> svDisruptions = germlineSV.getDisruptions();

            if(svDisruptions.isEmpty())
                return;

            for(Map.Entry<GeneAnnotation,String> entry : svDisruptions.entrySet())
            {
                final GeneAnnotation gene = entry.getKey();
                final String disruptionType = entry.getValue();
                final Transcript transcript = gene.canonical();

                mWriter.write(String.format("%s,%s,%s,%s,%s,%d,%d",
                        sampleId, germlineSV.Id, gene.isStart(), germlineSV.Type,
                        gene.chromosome(), gene.position(), gene.orientation()));

                mWriter.write(String.format(",%s,%s,%s,%d,%d,%s,%s,%s",
                        gene.StableId, gene.GeneName, transcript.StableId,
                        transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(), transcript.regionType(),
                        disruptionType));

                mWriter.newLine();
           }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing disruptions: {}", e.toString());
        }
    }

     */


    /* from main class - re-reading in SVs to annotate with disruptions

    private void reprocessVariants(final String sampleId, final List<GermlineSV> germlineSVs)
    {
        /*
        GM_LOGGER.info("sample({}) reprocessing {} germline SVs for disruptions", sampleId, germlineSVs.size());

        int svIndex = 0;
        for(final GermlineSV germlineSV : germlineSVs)
        {
            germlineSV.createSV();
            mGeneImpact.populateGeneAnnotations(germlineSV.sv(), svIndex++, germlineSV.getBreakendGenes(), germlineSV.getOverlappedGenes());
        }

        mGeneImpact.findDisruptiveVariants(sampleId, germlineSVs);

        int disruptiveSvCount = (int)germlineSVs.stream().filter(x -> !x.getDisruptions().isEmpty()).count();

        if(disruptiveSvCount > 0)
        {
            GM_LOGGER.debug("sample({}) has {} disruptive SVs from total({})",
                    germlineSVs.get(0).SampleId, disruptiveSvCount, germlineSVs.size());
        }

        germlineSVs.forEach(x -> writeCsv(x));
    }
    */
}
