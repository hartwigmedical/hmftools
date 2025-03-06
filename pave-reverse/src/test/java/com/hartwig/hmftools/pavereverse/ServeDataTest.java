package com.hartwig.hmftools.pavereverse;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Before;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ServeDataTest
{
    private final Map<String, Set<ServeItem>> geneToItems = new HashMap<>();
    //    private final File ensemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
    private final File ensemblDataDir = new File("/Users/timlavers/work/scratch/ensembl");
    private ReversePave baseSequenceVariantsCalculator;
    private Set<GeneAnnotation> knownDifferences;

    @Before
    public void setup() throws IOException
    {
        loadServeItems();
        createTransval();
        loadKnownDifferences();
    }

    private void createTransval() throws FileNotFoundException
    {
        final String genomePath = "/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";
        RefGenomeInterface genome = new RefGenomeSource(new IndexedFastaSequenceFile(new File(genomePath)));
        baseSequenceVariantsCalculator = new ReversePave(ensemblDataDir, RefGenomeVersion.V38, genome);
    }

    private void loadServeItems() throws IOException
    {
        File serveFile = new File("/Users/timlavers/work/junk/serve.json");
        String contents = Files.readString(serveFile.toPath());
        Gson gson = new Gson();
        var map = gson.fromJson(contents, JsonObject.class);
        var records = map.get("records");
        var hotspotsArray = ((JsonObject) records).get("V38")
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

    private void loadKnownDifferences() throws IOException
    {
        knownDifferences = new HashSet<>();
        knownDifferences.add(new GeneAnnotation("KDM4C", "G302E"));
        knownDifferences.add(new GeneAnnotation("KDM4C", "C778F"));
        knownDifferences.add(new GeneAnnotation("KDM4C", "K53M"));

        ClassLoader classLoader = getClass().getClassLoader();
        File file = new File(Objects.requireNonNull(classLoader.getResource("serve_differences.txt")).getFile());
        Files.readAllLines(file.toPath()).forEach(line ->
                {
                    String[] parts = line.split(" ");
                    knownDifferences.add(new GeneAnnotation(parts[0], parts[1]));
                }
        );
    }

    private static String str(JsonObject jsonObject, String key)
    {
        return jsonObject.get(key).getAsString();
    }

//    @Test
    public void check()
    {
        List<StatsForGene> results = new ArrayList<>();
        for(String gene : geneToItems.keySet())
        {
            results.add(checkItemsForGene(gene));
        }
        reportResults(results);
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

//    @Test
    public void examples()
    {
        BaseSequenceVariants variant = baseSequenceVariantsCalculator.calculateVariantAllowMultipleNonCanonicalTranscriptMatches("TP53", "Q331Q");
        Assert.assertEquals(6, variant.changes().size());
    }

    private StatsForGene checkItemsForGene(String gene)
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
        final StatsForGene statsForGene = new StatsForGene(baseSequenceVariantsCalculator, gene);
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
                    baseSequenceVariantsCalculator.calculateVariantAllowMultipleNonCanonicalTranscriptMatches(collator.mGene, collator.mAnnotation);
            return new VariantStatus(collator, baseSequenceVariants);
        }
        catch(Throwable t)
        {
            return VariantStatus.withProcessingException(collator, t);
        }
    }

    private void p(String msg)
    {
        System.out.println(msg);
    }
}

class VariantStatus
{
    static VariantStatus withProcessingException(@NotNull final ProteinAnnotationCollator collator, final @NotNull Throwable e)
    {
        VariantStatus status = new VariantStatus((collator));
        status.mProcessingError = e;
        return status;
    }

    @NotNull
    final ProteinAnnotationCollator collator;

    Exception mParseException = null;

    Throwable mProcessingError = null;

    BaseSequenceVariants variant;

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, BaseSequenceVariants variant)
    {
        this.collator = collator;
        this.variant = variant;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator)
    {
        this.collator = collator;
    }

    boolean usesNonCanonicalTranscript()
    {
        return !variant.mTranscript.IsCanonical;
    }

    boolean parsedOk()
    {
        return mParseException == null;
    }

    boolean hasProcessingError()
    {
        return mProcessingError != null;
    }

    boolean hotspotsSame()
    {
        var common = Sets.intersection(collator.hotspots, variant.changes());
        if(!common.isEmpty())
        {
            return true;
        }
        if(collator.mAnnotation.contains("ins"))
        {
            if(collator.hotspots.size() == 1 && variant.mChanges.size() == 1)
            {
                return hotspotsSameModuloAlt();
            }
        }

        if(collator.mAnnotation.endsWith("dup"))
        {
            if(collator.hotspots.size() == 1 && variant.mChanges.size() == 1)
            {
                return hotspotsSameModuloAlt();
            }

            return collator.hotspots.containsAll(variant.changes());
        }
        if(collator.mAnnotation.endsWith("fs"))
        {
            if(collator.hotspots.containsAll(variant.changes()))
            {
                return true;
            }
        }
        return false;
    }

    boolean hotspotsSameModuloAlt()
    {
        if(variant.changes().size() != 1)
        {
            return false;
        }
        if(collator.hotspots.size() != 1)
        {
            return false;
        }
        BaseSequenceChange calculatedHS = variant.changes().iterator().next();
        BaseSequenceChange givenHS = collator.hotspots.iterator().next();
        if(calculatedHS.mPosition != givenHS.mPosition)
        {
            return false;
        }
        if(!calculatedHS.Ref.equals(givenHS.Ref))
        {
            return false;
        }
        return givenHS.mChromosome.equals(calculatedHS.mChromosome);
    }

    boolean hotspotsSameModuloPosition()
    {
        Set<FloatingHotspot> variantFloatingHotspots = variant.changes().stream().map(FloatingHotspot::new).collect(Collectors.toSet());
        Set<FloatingHotspot> collatorFloatingHotspots = collator.hotspots.stream().map(FloatingHotspot::new).collect(Collectors.toSet());
        return variantFloatingHotspots.equals(collatorFloatingHotspots);
    }

    private BaseSequenceChange reduceHotspot(BaseSequenceChange hotspot)
    {
        DeletionInsertionChange delta = new DeletionInsertionChange(hotspot.Ref, hotspot.Alt);
        return delta.toHotspot(new ChangeLocation(hotspot.mChromosome, hotspot.mPosition));
    }
}

class FloatingHotspot
{
    @NotNull
    public final String mChromosome;

    @NotNull
    public final String Ref;

    @NotNull
    public final String Alt;

    public FloatingHotspot(BaseSequenceChange hotspot)
    {
        mChromosome = hotspot.mChromosome;
        Ref = hotspot.Ref;
        Alt = hotspot.Alt;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final FloatingHotspot that = (FloatingHotspot) o;
        return Objects.equals(mChromosome, that.mChromosome) && Objects.equals(Ref, that.Ref)
                && Objects.equals(Alt, that.Alt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Ref, Alt);
    }
}

class ServeItem
{
    @NotNull
    public final String Gene;

    @NotNull
    public final String Annotation;

    @NotNull
    public final String Chromosome;

    @NotNull
    public final String Ref;

    @NotNull
    public final String Alt;

    public final int Position;

    ServeItem(@NotNull final String gene, @NotNull final String annotation, @NotNull final String chromosome, @NotNull final String ref,
            @NotNull final String alt, final int position)
    {
        Preconditions.checkArgument(!gene.isEmpty());
        Preconditions.checkArgument(!annotation.isEmpty());
        Gene = gene;
        Annotation = annotation;
        Chromosome = chromosome;
        Ref = ref;
        Alt = alt;
        Position = position;
    }

    @Override
    public String toString()
    {
        return "ServeItem{" +
                "Gene='" + Gene + '\'' +
                ", Annotation='" + Annotation + '\'' +
                ", Chromosome='" + Chromosome + '\'' +
                ", Ref='" + Ref + '\'' +
                ", Alt='" + Alt + '\'' +
                ", Position=" + Position +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ServeItem serveItem = (ServeItem) o;
        return Position == serveItem.Position && Objects.equals(Gene, serveItem.Gene)
                && Objects.equals(Annotation, serveItem.Annotation) && Objects.equals(Chromosome, serveItem.Chromosome)
                && Objects.equals(Ref, serveItem.Ref) && Objects.equals(Alt, serveItem.Alt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Gene, Annotation, Chromosome, Ref, Alt, Position);
    }
}

class DifferenceWithTransvar {
    public final BaseSequenceVariants variant;
    public final ProteinAnnotationCollator collator;
    public final String type;

    DifferenceWithTransvar(final VariantStatus variantStatus, ReversePave baseSequenceVariantsCalculator)
    {
        this.variant = variantStatus.variant;
        this.collator = variantStatus.collator;
        ProteinVariant transvalVariant = baseSequenceVariantsCalculator.variationParser().parseGeneVariants(collator.mGene, collator.mAnnotation).iterator().next();
        type = transvalVariant.getClass().getSimpleName();
            System.out.println(type + " difference for: " + collator.mGene + " " + collator.mAnnotation
                    + ", calculated: " + variant.transcriptName()  + " canonical: " + variant.mTranscript.IsCanonical
                    + " " + variant.changes()
                    + " from transvar: " + collator.hotspots);
    }
}
class StatsForGene
{
    final ReversePave baseSequenceVariantsCalculator;
    @NotNull
    public final String mGene;

    public int NumberProcessedWithoutError = 0;
    public int NumberNotParsed = 0;

    public int NumberWithDifferentHotspots = 0;

    public int NumberWithSameHotspots = 0;

    public int NumberWithSameHotspotsThatUseNonCanonicalTranscript = 0;
    public int NumberWithDifferentHotspotsThatUseNonCanonicalTranscript = 0;
    public int NumberWithNoMatchingTranscript = 0;
    public int NumberForWhichNoVariantCouldBeCalculated = 0;
    public int NumberWithNoUniqueMatchingTranscript = 0;
    public int NumberWithProcessingError = 0;
    public int NumberWithChangesAcrossMultipleExons = 0;

    public final List<VariantStatus> AnnotationsWithDifferentHotspots = new ArrayList<>();
    public final List<DifferenceWithTransvar> Differences = new ArrayList<>();

    StatsForGene(final ReversePave baseSequenceVariantsCalculator, @NotNull final String gene)
    {
        this.baseSequenceVariantsCalculator = baseSequenceVariantsCalculator;
        mGene = gene;
    }

    public void recordResultsForAnnotation(VariantStatus comparison)
    {
        if(comparison.mParseException != null)
        {
            System.out.println(comparison.mParseException.getMessage());
        }
        if(comparison.hasProcessingError())
        {
            Throwable t = comparison.mProcessingError;
            System.out.println(t.getMessage());
        }
        if(comparison.hotspotsSame())
        {
            NumberWithSameHotspots++;
            if(comparison.usesNonCanonicalTranscript())
            {
                NumberWithSameHotspotsThatUseNonCanonicalTranscript++;
            }
        }
        else
        {
            AnnotationsWithDifferentHotspots.add(comparison);
            NumberWithDifferentHotspots++;
            Differences.add(new DifferenceWithTransvar(comparison, baseSequenceVariantsCalculator));
            if(comparison.usesNonCanonicalTranscript())
            {
                NumberWithDifferentHotspotsThatUseNonCanonicalTranscript++;
            }
        }
    }

    public void recordNotParsed()
    {
        NumberNotParsed++;
    }

    @Override
    public String toString()
    {
        return "StatsForGene{" +
                "Gene='" + mGene + '\'' +
                ", NumberNotParsed=" + NumberNotParsed +
                ", NumberWithDifferentHotspots=" + NumberWithDifferentHotspots +
                ", NumberWithSameHotspots=" + NumberWithSameHotspots +
                '}';
    }
}

class ProteinAnnotationCollator
{
    @NotNull
    final String mChromosome;

    @NotNull
    final String mGene;

    @NotNull
    final String mAnnotation;

    @NotNull
    Set<BaseSequenceChange> hotspots = new HashSet<>();

    public ProteinAnnotationCollator(ServeItem serveItem)
    {
        mChromosome = serveItem.Chromosome;
        mGene = serveItem.Gene;
        mAnnotation = serveItem.Annotation;
    }

    public void addHotspot(ServeItem serveItem)
    {
        Preconditions.checkArgument(serveItem.Chromosome.equals(mChromosome));
        Preconditions.checkArgument(serveItem.Gene.equals(mGene));
        Preconditions.checkArgument(serveItem.Annotation.equals(mAnnotation));
        hotspots.add(new BaseSequenceChange(serveItem.Ref, serveItem.Alt, mChromosome, serveItem.Position));
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ProteinAnnotationCollator that = (ProteinAnnotationCollator) o;
        return Objects.equals(mChromosome, that.mChromosome) && Objects.equals(mGene, that.mGene)
                && Objects.equals(mAnnotation, that.mAnnotation);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, mGene, mAnnotation);
    }

    @Override
    public String toString()
    {
        return "ProteinAnnotationCollator{" +
                "mChromosome='" + mChromosome + '\'' +
                ", mGene='" + mGene + '\'' +
                ", mAnnotation='" + mAnnotation + '\'' +
                '}';
    }
}

class GeneAnnotation implements Comparable<GeneAnnotation>
{
    @NotNull
    public final String mGene;

    @NotNull
    public final String mAnnotation;

    GeneAnnotation(@NotNull final String mGene, @NotNull final String mAnnotation)
    {
        this.mGene = mGene;
        this.mAnnotation = mAnnotation;
    }

    @Override
    public int compareTo(@NotNull final GeneAnnotation o)
    {
        int byGene = this.mGene.compareTo(o.mGene);
        if(byGene != 0)
        {
            return byGene;
        }
        return this.mAnnotation.compareTo(o.mAnnotation);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final GeneAnnotation that = (GeneAnnotation) o;
        return Objects.equals(mGene, that.mGene) && Objects.equals(mAnnotation, that.mAnnotation);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mGene, mAnnotation);
    }

    @Override
    public String toString()
    {
        return "GeneAnnotation{" +
                "mGene='" + mGene + '\'' +
                ", mAnnotation='" + mAnnotation + '\'' +
                '}';
    }
}

