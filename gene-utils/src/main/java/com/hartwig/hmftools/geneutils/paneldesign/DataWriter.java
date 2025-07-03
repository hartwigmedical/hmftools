package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NaN;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.tuple.Pair;

public class DataWriter
{
    public static final String PANEL_DEFINITION_FILE_EXTENSION = ".panel_definition.tsv";
    public static final String GENE_REGION_FILE_EXTENSION = ".gene_region.tsv";
    public static final String CANDIDATE_FILE_EXTENSION = ".probe_candidate.tsv";

    enum PanelDefinitionColumn
    {
        Chromosome,
        PositionStart,
        PositionEnd,
        RegionType,
        SourceInfo;
    }

    public static void writePanelDefinition(final String filename, final List<PanelRegion> panelRegions)
    {
        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            DelimFileWriter.write(writer, PanelDefinitionColumn.values(), panelRegions,
                    (r, row) -> {
                        row.set(PanelDefinitionColumn.Chromosome, r.Chromosome);
                        row.set(PanelDefinitionColumn.PositionStart, r.start());
                        row.set(PanelDefinitionColumn.PositionEnd, r.end());
                        row.set(PanelDefinitionColumn.RegionType, r.Type.toString());
                        row.set(PanelDefinitionColumn.SourceInfo, r.extendedSourceInfo());
             });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    enum GeneRegionColumn
    {
        GeneName,
        RegionType,
        Chromosome,
        RegionStart,
        RegionEnd,
        UseWholeRegion,
        ProbeStart,
        ProbeEnd,
        ProbeGcContent,
        ProbeQualityScore,
        ProbeSequence
    }

    public static void writeTargertedGeneRegions(final String filename, final List<TargetedGeneRegion> targetedGeneRegions)
    {
        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            DelimFileWriter.write(writer, GeneRegionColumn.values(), targetedGeneRegions,
                (r, row) -> {
                    row.set(GeneRegionColumn.GeneName, r.getGene().getGeneData().GeneName);
                    row.set(GeneRegionColumn.RegionType, r.getType().name());
                    row.set(GeneRegionColumn.Chromosome, r.getChromosome());
                    row.set(GeneRegionColumn.RegionStart, r.getStart());
                    row.set(GeneRegionColumn.RegionEnd, r.getEnd());
                    row.set(GeneRegionColumn.UseWholeRegion, r.useWholeRegion());

                    // some we use whole region
                    if(!r.useWholeRegion())
                    {
                        ProbeCandidate selectedProbe = r.getSelectedProbe();
                        if(selectedProbe != null)
                        {
                            row.set(GeneRegionColumn.ProbeStart, selectedProbe.getStart());
                            row.set(GeneRegionColumn.ProbeEnd, selectedProbe.getEnd());
                            row.set(GeneRegionColumn.ProbeGcContent, selectedProbe.getGcContent());
                            row.set(GeneRegionColumn.ProbeQualityScore, selectedProbe.getQualityScore().orElse(NaN));
                        }
                    }
                });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    enum GeneProbeCandidateColumn
    {
        GeneName,
        RegionType,
        Chromosome,
        RegionStart,
        RegionEnd,
        ProbeStart,
        ProbeEnd,
        ProbeGcContent,
        ProbeQualityScore,
        Selected,
        ProbeSequence
    }

    public static void writeCandidates(final String filename, final List<TargetedGeneRegion> targetedGeneRegions)
    {
        List<Pair<TargetedGeneRegion, ProbeCandidate>> probeList = targetedGeneRegions.stream()
                .flatMap(o -> o.getProbeCandidates().stream().map(c -> Pair.of(o, c)))
                .collect(Collectors.toList());

        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            DelimFileWriter.write(writer, GeneProbeCandidateColumn.values(), probeList,
                    (r, row) -> {
                        row.set(GeneProbeCandidateColumn.GeneName, r.getLeft().getGene().getGeneData().GeneName);
                        row.set(GeneProbeCandidateColumn.RegionType, r.getLeft().getType().name());
                        row.set(GeneProbeCandidateColumn.Chromosome, r.getLeft().getChromosome());
                        row.set(GeneProbeCandidateColumn.RegionStart, r.getLeft().getStart());
                        row.set(GeneProbeCandidateColumn.RegionEnd, r.getLeft().getEnd());
                        row.set(GeneProbeCandidateColumn.ProbeStart, r.getRight().getStart());
                        row.set(GeneProbeCandidateColumn.ProbeEnd, r.getRight().getEnd());
                        row.set(GeneProbeCandidateColumn.ProbeGcContent, r.getRight().getGcContent());
                        row.set(GeneProbeCandidateColumn.ProbeQualityScore, r.getRight().getQualityScore().orElse(NaN));
                        row.set(GeneProbeCandidateColumn.Selected, Boolean.toString(r.getLeft().getSelectedProbe() == r.getRight()));
                        row.set(GeneProbeCandidateColumn.ProbeSequence, r.getRight().getSequence());
                    });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

}
