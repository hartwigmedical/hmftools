package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

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

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DisruptionFinder
{
    private Set<String> mDisruptionGeneIDPanel;
    private BufferedWriter mWriter;
    private final String mOutputDir;

    public static final String DRUP_TSG_GENES_FILE = "drup_tsg_file";

    private static final Logger LOGGER = LogManager.getLogger(DisruptionFinder.class);

    public DisruptionFinder(final CommandLine cmd, final SvGeneTranscriptCollection geneTransCache, final String outputDir)
    {
        mDisruptionGeneIDPanel = null;
        mOutputDir = outputDir;
        mWriter = null;

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

    public void writeSampleData(final String sampleId, final List<Transcript> disruptions)
    {
        // write sample file for patient reporter
        List<ReportableDisruption> reportedDisruptions = Lists.newArrayList();

        for(final Transcript transcript : disruptions)
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

    public void writeDisruptionData(final String sampleId, final Transcript transcript)
    {
        if(mWriter == null)
            return;

        try
        {
            final GeneAnnotation gene = transcript.parent();

            mWriter.write(String.format("%s,%s,%d,%s,%s,%d,%d",
                    sampleId, transcript.reportableDisruption(), gene.id(), gene.isStart(),
                    gene.chromosome(), gene.position(), gene.orientation()));

            mWriter.write(String.format(",%s,%s,%d,%s,%d,%d,%s,%s",
                    gene.StableId, gene.GeneName, gene.Strand, transcript.StableId,
                    transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(), transcript.regionType()));

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
