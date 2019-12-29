package com.hartwig.hmftools.linx.vcf_filters;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.GENE_PANEL_FILE;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.GENE_TRANSCRIPTS_DIR;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneImpact
{
    private final GermlineVcfConfig mConfig;
    private final SvGeneTranscriptCollection mGeneCollection;
    private final Map<String, List<EnsemblGeneData>> mGenePanel;

    // map of SV Id to a pair of gene annotations (for start and end breakends)
    private final Map<String,List<List<GeneAnnotation>>> mSvGeneAnnotations;

    private final Map<StructuralVariant,List<GeneAnnotation>> mSvGeneDisruptions;

    private static final Logger LOGGER = LogManager.getLogger(GeneImpact.class);

    public GeneImpact(final GermlineVcfConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mGenePanel = Maps.newHashMap();
        mSvGeneAnnotations = Maps.newHashMap();
        mSvGeneDisruptions = Maps.newHashMap();

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            mGeneCollection = new SvGeneTranscriptCollection();
            mGeneCollection.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));
            mGeneCollection.loadEnsemblData(true);

            mGeneCollection.setRequiredData(config.CheckDisruptions, false, false, true);

            if(config.CheckDisruptions)
            {
                mGeneCollection.loadEnsemblTranscriptData(Lists.newArrayList());
            }

            if(cmd.hasOption(GENE_PANEL_FILE))
                loadGenePanel(cmd.getOptionValue(GENE_PANEL_FILE));
        }
        else
        {
            mGeneCollection = null;
        }
    }

    public final Map<StructuralVariant,List<GeneAnnotation>> getGeneDisruptions() { return mSvGeneDisruptions; }

    private void loadGenePanel(final String geneFile)
    {
        // gene_panel.csv
        if (!Files.exists(Paths.get(geneFile)))
            return;

        try
        {
            List<String> fileContents = Files.readAllLines(new File(geneFile).toPath());

            for(final String gene : fileContents)
            {
                EnsemblGeneData geneData = mGeneCollection.getGeneDataByName(gene);

                if(geneData == null)
                {
                    LOGGER.warn("gene({}) not found in Ensembl data cache", gene);
                    continue;
                }

                List<EnsemblGeneData> chrGenesList = mGenePanel.get(geneData.Chromosome);

                if(chrGenesList == null)
                {
                    chrGenesList = Lists.newArrayList();
                    mGenePanel.put(geneData.Chromosome, chrGenesList);
                }

                chrGenesList.add(geneData);
            }

            LOGGER.info("loaded genePanelFile({}) with {} genes", geneFile, fileContents.size());
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load gene panel file({}): {}", geneFile, e.toString());
        }
    }

    public void annotationWithGenes(final StructuralVariant sv, final String[] geneAnnotations)
    {
        if(mGeneCollection == null)
            return;

        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(sv.type() == SGL && be == SE_END)
                continue;

            final List<EnsemblGeneData> matchedGenes = mGeneCollection.findGenes(sv.chromosome(isStart(be)), sv.position(isStart(be)), 0);

            if(!matchedGenes.isEmpty())
            {
                String genesStr = "";

                for(final EnsemblGeneData gene : matchedGenes)
                {
                    genesStr = appendStr(genesStr, gene.GeneName, ';');
                }

                geneAnnotations[be] = genesStr;
            }
        }
    }

    public String annotateWithGenePanel(final StructuralVariant sv)
    {
        if(sv.type() != DEL && sv.type() != DUP)
            return "";

        List<EnsemblGeneData> genesList = mGenePanel.get(sv.chromosome(true));

        if(genesList == null)
            return "";

        String genesStr = "";

        for(final EnsemblGeneData geneData : genesList)
        {
            if(sv.position(true) < geneData.GeneStart && sv.position(false) > geneData.GeneEnd)
            {
                genesStr = appendStr(genesStr, geneData.GeneName, ';');
            }
        }

        return genesStr;
    }

    public boolean hasGeneAnnotation(final StructuralVariant sv)
    {
        List<EnsemblGeneData> genesList = mGenePanel.get(sv.chromosome(true));

        if(genesList != null)
        {
            for (final EnsemblGeneData geneData : genesList)
            {
                // fully overlapping DEL or DUP
                if (sv.type() == DEL || sv.type() == DUP)
                {
                    if (sv.position(true) < geneData.GeneStart && sv.position(false) > geneData.GeneEnd)
                        return true;
                }

                // breakend falling in the gene
                if (sv.position(true) > geneData.GeneStart && sv.position(true) < geneData.GeneEnd)
                    return true;

                if (sv.type() == DEL || sv.type() == DUP | sv.type() == INV)
                {
                    if (sv.position(false) > geneData.GeneStart && sv.position(false) < geneData.GeneEnd)
                        return true;
                }
            }
        }

        // check other end of BND
        if(sv.type() == BND)
        {
            genesList = mGenePanel.get(sv.chromosome(false));

            if(genesList == null)
                return false;

            for (final EnsemblGeneData geneData : genesList)
            {
                if (sv.position(false) > geneData.GeneStart && sv.position(false) < geneData.GeneEnd)
                    return true;
            }
        }

        return false;
    }

    public void findDisruptiveVariants(final List<StructuralVariant> svList)
    {
        // SVs are not disruptive if they:
        // - are outside a gene
        // - form a TI within a single intron and the other breakend is outside a gene
        // - are a DUP or DEL and contained within an intron

        mSvGeneDisruptions.clear();

        for(int index = 0; index < svList.size(); ++index)
        {
            final StructuralVariant var = svList.get(index);
            List<List<GeneAnnotation>> breakendPairGenesList = populateGeneAnnotations(var, index);

            if(breakendPairGenesList == null)
                continue;

            final List<GeneAnnotation> genesStart = breakendPairGenesList.get(SE_START);
            final List<GeneAnnotation> genesEnd = breakendPairGenesList.get(SE_END);

            if(genesStart.isEmpty() || genesEnd.isEmpty())
            {
                final List<GeneAnnotation> genesList = !genesStart.isEmpty() ? genesStart : genesEnd;

                addGeneDisruption(var, genesList);

                // one breakend disrupts a gene
                final Transcript trans = genesList.get(0).canonical();

                LOGGER.debug("var({}:{} {}:{} -> {}:{}) gene({}) disrupted by {} breakend exons({} & {})",
                        var.id(), var.type(), var.chromosome(true), var.position(true),
                        var.end() != null ? var.chromosome(false) : "0", var.end() != null ? var.position(false) : -1,
                        trans.geneName(), !genesStart.isEmpty() ? "start" : "end", trans.ExonUpstream, trans.ExonDownstream);
            }
            else if(var.type() == BND || var.type() == INV)
            {
                // both ends are in genes and initially are considered disruptive
                addGeneDisruption(var, genesStart);
                addGeneDisruption(var, genesEnd);

                final Transcript transStart = !genesStart.isEmpty() ? genesStart.get(0).canonical() : null;
                final Transcript transEnd = !genesEnd.isEmpty() ? genesEnd.get(0).canonical() : null;

                LOGGER.debug("var({}:{} {}:{} -> {}:{}) genes({} & {}) disrupted",
                        var.id(), var.type(), var.chromosome(true), var.position(true),
                        var.end() != null ? var.chromosome(false) : "0", var.end() != null ? var.position(false) : -1,
                        transStart != null ? transStart.geneName() : "none", transEnd != null ? transEnd.geneName() : "none");
            }
            else if(isDisruptiveSV(var, genesStart, genesEnd))
            {
                addGeneDisruption(var, genesStart);
                addGeneDisruption(var, genesEnd);
                continue;
            }
        }

        // check for SVs linked within the same intron
        removeIntronicLinks();
    }

    private void addGeneDisruption(final StructuralVariant var, final List<GeneAnnotation> genesList)
    {
        List<GeneAnnotation> geneDisruptions = mSvGeneDisruptions.get(var);

        if(geneDisruptions == null)
        {
            geneDisruptions = Lists.newArrayList();
            mSvGeneDisruptions.put(var, geneDisruptions);
        }

        geneDisruptions.addAll(genesList);
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

    private void removeIntronicLinks()
    {
        final List<StructuralVariant> intronicLinkedSVs = Lists.newArrayList();

        for(Map.Entry<StructuralVariant,List<GeneAnnotation>> entry1 : mSvGeneDisruptions.entrySet())
        {
            final StructuralVariant var1 = entry1.getKey();

            if(intronicLinkedSVs.contains(var1))
                continue;

            final List<GeneAnnotation> genesList1 = entry1.getValue();

            // ignore SVs with both breakends in genes since any possibility of chaining is too difficult to establish
            boolean dualBreakends1 = genesList1.stream().anyMatch(x -> x.isStart()) && genesList1.stream().anyMatch(x -> !x.isStart());

            if(dualBreakends1 && var1.type() != INV)
                continue;

            for (Map.Entry<StructuralVariant, List<GeneAnnotation>> entry2 : mSvGeneDisruptions.entrySet())
            {
                final StructuralVariant var2 = entry2.getKey();

                if(intronicLinkedSVs.contains(var2))
                    continue;

                if(var1.equals(var2))
                    continue;

                final List<GeneAnnotation> genesList2 = entry2.getValue();

                boolean dualBreakends2 = genesList2.stream().anyMatch(x -> x.isStart()) && genesList2.stream().anyMatch(x -> !x.isStart());

                if(dualBreakends2 && var2.type() != INV)
                    continue;

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
                            long tiLength = abs(gene2.position() - gene1.position());

                            LOGGER.info("SVs({} & {}) have facing intronic breakends in gene({}) exons({} -> {}) tiLength({})",
                                    var1.id(), var2.id(), gene1.GeneName, trans1.ExonUpstream, trans1.ExonDownstream, tiLength);

                            intronicLinkedSVs.add(var1);
                            intronicLinkedSVs.add(var2);
                        }
                    }
                }
                else if(dualBreakends1 && dualBreakends2)
                {
                    final Transcript transStart1 = genesList1.stream().filter(x -> x.isStart()).findFirst().orElse(null).canonical();
                    final Transcript transEnd1 = genesList1.stream().filter(x -> !x.isStart()).findFirst().orElse(null).canonical();
                    final Transcript transStart2 = genesList2.stream().filter(x -> x.isStart()).findFirst().orElse(null).canonical();
                    final Transcript transEnd2 = genesList2.stream().filter(x -> !x.isStart()).findFirst().orElse(null).canonical();

                    if(!transStart1.isIntronic() || !transEnd1.isIntronic() || !transStart2.isIntronic() || !transEnd2.isIntronic())
                        continue;

                    // check orientations are opposite
                    if(transStart1.gene().orientation() == transStart2.gene().orientation())
                        continue;

                    if(transStart1.geneName().equals(transEnd1.geneName()) && transStart1.ExonDownstream == transEnd1.ExonDownstream
                    && transStart1.geneName().equals(transStart2.geneName()) && transStart1.ExonDownstream == transStart2.ExonDownstream
                    && transStart1.geneName().equals(transEnd2.geneName()) && transStart1.ExonDownstream == transEnd2.ExonDownstream)
                    {
                        long dbLength = abs(transStart1.gene().position() - transStart2.gene().position());

                        LOGGER.info("SVs({} & {}) form reciprocal INV in gene({}) exons({} -> {}) dbLength({})",
                                var1.id(), var2.id(), transStart1.geneName(),
                                transStart1.ExonUpstream, transStart1.ExonDownstream, dbLength);

                        intronicLinkedSVs.add(var1);
                        intronicLinkedSVs.add(var2);
                    }
                }
            }
        }

        intronicLinkedSVs.stream().forEach(x -> mSvGeneDisruptions.remove(x));
    }

    private Transcript findMatchingTrans(
            final List<StructuralVariant> svList, final List<Integer> disruptiveSVs,
            final Transcript otherTrans, final StructuralVariant otherrVar)
    {

        return null;
    }


    private List<List<GeneAnnotation>> populateGeneAnnotations(final StructuralVariant var, int svIndex)
    {
        List<List<GeneAnnotation>> breakendPairGenesList = null;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            boolean isStart = isStart(se);

            List<GeneAnnotation> genesList = mGeneCollection.findGeneAnnotationsBySv(
                    svIndex, isStart, var.chromosome(isStart), var.position(isStart), var.orientation(isStart), 0);

            // check for presence of a canonical transcript for each gene
            genesList = genesList.stream().filter(x -> x.canonical() != null).collect(Collectors.toList());

            if(genesList.isEmpty())
                continue;

            genesList.stream().forEach(x -> x.setPositionalData(var.chromosome(isStart), var.position(isStart), var.orientation(isStart)));

            if(breakendPairGenesList == null)
            {
                breakendPairGenesList = Lists.newArrayList();
                breakendPairGenesList.add(Lists.newArrayList());
                breakendPairGenesList.add(Lists.newArrayList());
                mSvGeneAnnotations.put(var.id(), breakendPairGenesList);
            }


            List<GeneAnnotation> svGenes = breakendPairGenesList.get(se);
            svGenes.addAll(genesList);
        }

        return breakendPairGenesList;
    }

}
