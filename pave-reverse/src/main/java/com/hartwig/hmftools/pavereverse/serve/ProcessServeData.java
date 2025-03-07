package com.hartwig.hmftools.pavereverse.serve;

import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.RPV_LOGGER;
import static com.hartwig.hmftools.pavereverse.batch.VariantsEncoder.columns;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ReversePave;
import com.hartwig.hmftools.pavereverse.batch.VariantRow;
import com.hartwig.hmftools.pavereverse.batch.VariantsEncoder;

import org.apache.logging.log4j.util.BiConsumer;
import org.jetbrains.annotations.NotNull;

public class ProcessServeData
{
    private final Map<String, Set<ServeItem>> geneToItems = new HashMap<>();
    private final ReversePave reversePave;
    private final RefGenomeVersion mRefGenomeVersion;

    public ProcessServeData(final ReversePave reversePave, final RefGenomeVersion refGenomeVersion)
    {
        this.reversePave = reversePave;
        this.mRefGenomeVersion = refGenomeVersion;
    }

    public void checkServeData(String serveFilePath, String outputTsvFilePath) throws IOException
    {
        loadServeItems(serveFilePath);
        BiConsumer<VariantRow, DelimFileWriter.Row> encoder = new VariantsEncoder();
        DelimFileWriter<VariantRow> writer = new DelimFileWriter<>(outputTsvFilePath, columns(), encoder);
        List<StatsForGene> results = new ArrayList<>();
        for(String gene : geneToItems.keySet())
        {
            results.add(processVariantsForGene(gene, writer));
        }
        reportResults(results);
        writer.close();
    }

    private void loadServeItems(String serveFilePath) throws IOException
    {
        String version = "V" + mRefGenomeVersion.identifier();
        File serveFile = new File(serveFilePath);
        String contents = Files.readString(serveFile.toPath());
        Gson gson = new Gson();
        var map = gson.fromJson(contents, JsonObject.class);
        var records = map.get("records");
        var hotspotsArray = ((JsonObject) records).get(version)
                .getAsJsonObject()
                .get("knownEvents")
                .getAsJsonObject().get("hotspots").getAsJsonArray();
        int start = 0;
        int stop = hotspotsArray.size();
        for(int i = start; i < stop; i++)
        {
            JsonObject hotspot = hotspotsArray.get(i).getAsJsonObject();
            final String gene = str(hotspot, "gene");
            final String annotation = str(hotspot, "inputProteinAnnotation");
            final String chromosome = str(hotspot, "chromosome");
            final String ref = str(hotspot, "ref");
            final String alt = str(hotspot, "alt");
            final int position = hotspot.get("position").getAsInt();
            if(annotation.isBlank())
            {
                continue;
            }
            ServeItem item = new ServeItem(gene, annotation, chromosome, ref, alt, position);
            if(!geneToItems.containsKey(gene))
            {
                geneToItems.put(gene, new HashSet<>());
            }
            geneToItems.get(gene).add(item);
        }
    }

    private static String str(JsonObject jsonObject, String key)
    {
        return jsonObject.get(key).getAsString();
    }

    private void reportResults(@NotNull List<StatsForGene> statsForGenes)
    {
        int numberOfGenes = 0;
        int numberNotParsed = 0;
        int numberWithSameHotspots = 0;
        int numberWithNoTranscriptFound = 0;
        int numberWithMultipleTranscriptFound = 0;
        int numberWithProcessingError = 0;
        int numberWithVariantsAcrossMultipleExons = 0;
        int numberProcessedWithoutError = 0;
        int numberForWhichNoVariantCouldBeCalculated = 0;
        Map<String, Integer> differenceTypes = new HashMap<>();
        for(StatsForGene stats : statsForGenes)
        {
            numberOfGenes++;
            numberNotParsed += stats.NumberNotParsed;
            numberWithSameHotspots += stats.NumberWithSameHotspots;
            numberWithNoTranscriptFound += stats.NumberWithNoMatchingTranscript;
            numberWithMultipleTranscriptFound += stats.NumberWithNoUniqueMatchingTranscript;
            numberProcessedWithoutError += stats.NumberProcessedWithoutError;
            numberWithProcessingError += stats.NumberWithProcessingError;
            numberWithVariantsAcrossMultipleExons += stats.NumberWithChangesAcrossMultipleExons;
            numberForWhichNoVariantCouldBeCalculated += stats.NumberForWhichNoVariantCouldBeCalculated;
            stats.Differences.forEach(d -> {
                Integer typeCount = differenceTypes.remove(d.type);
                Integer newCount = typeCount == null ? 1 : typeCount + 1;
                differenceTypes.put(d.type, newCount);
            });
        }
        p("OVERALL STATS");
        p("*************");
        p("Number of genes: " + numberOfGenes);
        p("Number of annotations not parsed: " + numberNotParsed);
        p("Number with no transcript found: " + numberWithNoTranscriptFound);
        p("Number with multiple non-canonical transcripts giving different results: " + numberWithMultipleTranscriptFound);
        p("Number with variants across multiple exons (not handled): " + numberWithVariantsAcrossMultipleExons);
        p("Number processed without error: " + numberProcessedWithoutError);
        p("Number for which no variant could be calculate: " + numberForWhichNoVariantCouldBeCalculated);
        p("Number with processing error: " + numberWithProcessingError);
        p("Number of annotations with similar hotspots: " + numberWithSameHotspots);
        p("Difference types: " + differenceTypes);
    }

    private StatsForGene processVariantsForGene(String gene, DelimFileWriter<VariantRow> outputWriter)
    {
        Map<String, ProteinAnnotationCollator> collators = new HashMap<>();
        geneToItems.get(gene).forEach(serveItem ->
        {
            if(!collators.containsKey(serveItem.Annotation))
            {
                collators.put(serveItem.Annotation, new ProteinAnnotationCollator(serveItem));
            }
            collators.get(serveItem.Annotation).addHotspot(serveItem);
        });
        final StatsForGene statsForGene = new StatsForGene(reversePave, gene);
        collators.keySet().forEach(annotation ->
        {
            ProteinAnnotationCollator collator = collators.get(annotation);
            VariantStatus comparison = checkVariant(collator);
            if(comparison.parsedOk())
            {
                if(comparison.hasProcessingError())
                {
                    String msg = comparison.mProcessingError.getMessage();
                    if(msg == null)
                    {
                        p("This one had an error with no message: " + collator);
                        p("The processing error was: " + comparison.mProcessingError);
                        statsForGene.NumberWithProcessingError++;
                    }
                    else
                    {
                        if(msg.startsWith("No transcript found for gene"))
                        {
                            statsForGene.NumberWithNoMatchingTranscript++;
                        }
                        else if(msg.endsWith("but produce different hotspots"))
                        {
                            statsForGene.NumberWithNoUniqueMatchingTranscript++;
                        }
                        else if(msg.startsWith("Window end overlaps by more than one codon"))
                        {
                            statsForGene.NumberWithChangesAcrossMultipleExons++;
                        }
                        else if(msg.startsWith("No variant could be calculated"))
                        {
                            statsForGene.NumberForWhichNoVariantCouldBeCalculated++;
                        }
                        else
                        {
                            p("Processing error for: " + collator);
                            statsForGene.NumberWithProcessingError++;
                            p(msg);
                        }
                    }
                }
                else
                {
                    statsForGene.recordResultsForAnnotation(comparison);
                    statsForGene.NumberProcessedWithoutError++;
                    outputWriter.writeRow(new VariantRow(comparison.collator.mGene, comparison.collator.mAnnotation, comparison.variant));
                }
            }
            else
            {
                p(gene + " " + annotation);
                statsForGene.recordNotParsed();
            }
        });
        return statsForGene;
    }

    private VariantStatus checkVariant(ProteinAnnotationCollator collator)
    {
        try
        {
            BaseSequenceVariants baseSequenceVariants =
                    reversePave.calculateVariantAllowMultipleNonCanonicalTranscriptMatches(collator.mGene, collator.mAnnotation);
            return new VariantStatus(collator, baseSequenceVariants);
        }
        catch(Throwable t)
        {
            return VariantStatus.withProcessingException(collator, t);
        }
    }

    private void p(String msg)
    {
        RPV_LOGGER.info(msg);
    }
}
