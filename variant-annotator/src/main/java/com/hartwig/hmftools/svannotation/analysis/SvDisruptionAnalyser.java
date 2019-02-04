package com.hartwig.hmftools.svannotation.analysis;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableGeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SvDisruptionAnalyser
{
    private final Set<String> mDisruptionGeneIDPanel;
    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(SvDisruptionAnalyser.class);

    public SvDisruptionAnalyser()
    {
        mDisruptionGeneIDPanel = tsgDriverGeneIDs();
        mWriter = null;
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

    public final List<GeneDisruption> findDisruptions(final List<StructuralVariantAnnotation> annotations)
    {
        LOGGER.debug("finding disruptions in {} annotations", annotations.size());

        final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

        for (final StructuralVariantAnnotation annotation : annotations)
        {
            final boolean pureIntronicDisruptionCanonical = annotation.start()
                    .stream()
                    .filter(gene -> gene.canonical() != null)
                    .anyMatch(gene -> annotation.end()
                            .stream()
                            .filter(other -> other.canonical() != null)
                            .anyMatch(other -> intronicDisruptionOnSameTranscript(gene.canonical(), other.canonical())));

            if (pureIntronicDisruptionCanonical && annotation.variant().type() != StructuralVariantType.INV)
                continue;

            geneAnnotations.addAll(annotation.annotations());
        }

        final Multimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        geneAnnotations.forEach(g -> geneMap.put(g.GeneName, g));

        final List<GeneDisruption> disruptions = Lists.newArrayList();

        for (final String geneName : geneMap.keySet())
        {
            for (final GeneAnnotation gene : geneMap.get(geneName))
            {
                for (final Transcript transcript : gene.transcripts())
                {
                    final GeneDisruption disruption = ImmutableGeneDisruption.builder()
                            .reportable(mDisruptionGeneIDPanel.stream().anyMatch(geneID -> gene.synonyms().contains(geneID))
                                    && transcript.isCanonical())
                            .linkedAnnotation(transcript)
                            .build();

                    disruptions.add(disruption);
                }
            }
        }

        setBreakendDisruptive(annotations);

        return disruptions;
    }

    private static boolean intronicDisruptionOnSameTranscript(@NotNull Transcript t1, @NotNull Transcript t2)
    {
        boolean sameTranscript = t1.StableId.equals(t2.StableId);
        boolean bothIntronic = t1.isIntronic() && t2.isIntronic();
        boolean sameExonUpstream = t1.exonUpstream() == t2.exonUpstream();

        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    public static void setBreakendDisruptive(final List<StructuralVariantAnnotation> annotations)
    {
        // A breakend is not disruptive in a transcript if there exists another breakend where:
        // transcriptId = transcriptId AND svId = svId  AND
        //  isStartEnd <> isStartEnd AND ExonRankUpstream = ExonRankUpstream AND
        // sv.Type in (‘DEL’ ,’INS’, ‘DUP’)

        for(final StructuralVariantAnnotation annotation : annotations)
        {
            for(final GeneAnnotation startGene : annotation.start())
            {
                for (final Transcript trans1 : startGene.transcripts())
                {
                    for (final GeneAnnotation endGene : annotation.end())
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

        if(trans1.exonUpstream() != trans2.exonUpstream())
            return true;

        return false;
    }

    public void writeDisruptions(final List<GeneDisruption> disruptions, final String outputDir, final String sampleId, boolean hasMultipleSamples)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = outputDir;

                if (!outputFilename.endsWith(File.separator))
                    outputFilename += File.separator;

                if(hasMultipleSamples)
                    outputFilename += "DISRUPTIONS.csv";
                else
                    outputFilename += sampleId + "_disruptions.csv";

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,Reportable,SvId,Chromosome,Position,Orientation,Type,Ploidy");
                mWriter.write(",Gene,ChrBand,Transcript,Strand,RegionType,CodingType,Canonical,Biotype,ExonUp,ExonDown,IsDisruptive");
                mWriter.newLine();
            }

            BufferedWriter writer = mWriter;

            for(final GeneDisruption disruption : disruptions)
            {
                final Transcript transcript = disruption.linkedAnnotation();

                final GeneAnnotation svBreakend = transcript.parent();

                writer.write(String.format("%s,%s,%d,%s,%d,%d,%s,%.6f",
                        sampleId, disruption.reportable(), svBreakend.id(),
                        svBreakend.chromosome(), svBreakend.position(), svBreakend.orientation(),
                        svBreakend.type(), svBreakend.ploidy()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%s,%s,%s,%s,%d,%d,%s",
                                svBreakend.GeneName, svBreakend.karyotypeBand(), transcript.StableId, svBreakend.Strand,
                                transcript.regionType(), transcript.codingType(), transcript.isCanonical(), transcript.bioType(),
                                transcript.exonUpstream(), transcript.exonDownstream(), transcript.isDisruptive()));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing disruptions: {}", e.toString());
        }
    }

    public void onCompleted()
    {
        closeBufferedWriter(mWriter);
    }

}
