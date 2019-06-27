package com.hartwig.hmftools.linx.fusion;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DisruptionFinder
{
    private Set<String> mDisruptionGeneIDPanel;

    public static final String DRUP_TSG_GENES_FILE = "drup_tsg_file";

    private static final Logger LOGGER = LogManager.getLogger(DisruptionFinder.class);

    public DisruptionFinder(final CommandLine cmd, final SvGeneTranscriptCollection geneTransCache)
    {
        mDisruptionGeneIDPanel = null;

        initialise(cmd, geneTransCache);
    }

    private void initialise(final CommandLine cmd, final SvGeneTranscriptCollection geneTransCache)
    {
        mDisruptionGeneIDPanel = tsgDriverGeneIDs();

        // TEMP: load DRUP TSGs from file
        if(cmd != null && cmd.hasOption(DRUP_TSG_GENES_FILE))
        {
            loadDrupTSGs(cmd.getOptionValue(DRUP_TSG_GENES_FILE), geneTransCache);
        }
    }

    public boolean matchesDisruptionGene(final GeneAnnotation gene)
    {
        return mDisruptionGeneIDPanel.stream().anyMatch(geneID -> gene.synonyms().contains(geneID));
    }

    public static void markNonDisruptiveTranscripts(List<GeneAnnotation> genesStart, List<GeneAnnotation> genesEnd)
    {
        if(genesStart.isEmpty() || genesEnd.isEmpty())
            return;

        for(final GeneAnnotation startGene : genesStart)
        {
            for (final Transcript trans1 : startGene.transcripts())
            {
                for (final GeneAnnotation endGene : genesEnd)
                {
                    for (final Transcript trans2 : endGene.transcripts())
                    {
                        if(!areDisruptivePair(trans1, trans2))
                        {
                            trans1.setIsDisruptive(false);
                            trans2.setIsDisruptive(false);
                        }
                    }
                }
            }
        }
    }


    public static boolean areDisruptivePair(final Transcript trans1, final Transcript trans2)
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

    private void loadDrupTSGs(final String filename, final SvGeneTranscriptCollection geneTransCache)
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

                final EnsemblGeneData geneData = geneTransCache.getGeneDataByName(geneName);
                if(geneData != null)
                {
                    mDisruptionGeneIDPanel.add(geneData.GeneId);
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

}
