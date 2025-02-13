package com.hartwig.hmftools.pave.transval;

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
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class ServeDataTest
{
    private final Map<String, Set<ServeItem>> geneToItems = new HashMap<>();
    //    private final File ensemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
    private final File ensemblDataDir = new File("/Users/timlavers/work/scratch/ensembl");
    private Transval transval;
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
        transval = new Transval(ensemblDataDir, genome);
    }

    private void loadServeItems() throws IOException
    {
        File serveFile = new File("/Users/timlavers/work/junk/serve.json");
        String contents = Files.readString(serveFile.toPath());
        Gson gson = new Gson();
        var map = gson.fromJson(contents, JsonObject.class);
        var records = map.get("records");
        var hotspotsArray =
                ((JsonObject) records).get("V38").getAsJsonObject().get("knownEvents").getAsJsonObject().get("hotspots").getAsJsonArray();
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
//                        if (gene.equals("FGFR1") && annotation.equals("P283T"))
            ServeItem item = new ServeItem(gene, annotation, chromosome, ref, alt, position);
            if(!annotation.contains("ins")) continue;
            if(annotation.contains("del")) continue;
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
//        int numberWithSameHotspotsAndUsingNonCanonicalTranscript = 0;
//        int numberWithDifferentHotspotsAndUsingNonCanonicalTranscript = 0;
        SortedSet<GeneAnnotation> unaccountedDifferences = new TreeSet<>();
        SortedSet<GeneAnnotation> unaccountedDifferencesButWithHotspotsSameApartFromPosition = new TreeSet<>();
        for(StatsForGene stats : statsForGenes)
        {
            numberOfGenes++;
            numberNotParsed += stats.NumberNotParsed;
            numberWithSameHotspots += stats.NumberWithSameHotspots;
//            numberWithSameHotspotsAndUsingNonCanonicalTranscript += stats.NumberWithSameHotspotsThatUseNonCanonicalTranscript;
//            numberWithDifferentHotspotsAndUsingNonCanonicalTranscript += stats.NumberWithDifferentHotspotsThatUseNonCanonicalTranscript;
            stats.AnnotationsWithDifferentHotspots.forEach(variantStatus ->
            {
                if(!knownDifferences.contains(variantStatus.geneAnnotation()))
                {
                    unaccountedDifferences.add(variantStatus.geneAnnotation());
                    if (variantStatus.hotspotsSameModuloPosition())
                    {
                        unaccountedDifferencesButWithHotspotsSameApartFromPosition.add(variantStatus.geneAnnotation());
                    }
                }
            });
        }
        p("UNNACOUNTED DIFFERENCES:");
        p("************************");
        unaccountedDifferences.forEach(difference ->
                p(difference.mGene + " " + difference.mAnnotation));
        p("OVERALL STATS");
        p("*************");
        p("Number of genes: " + numberOfGenes);
        p("Number of annotations not parsed: " + numberNotParsed);
        p("Number of annotations with same hotspots: " + numberWithSameHotspots);
        p("Unaccounted differences: " + unaccountedDifferences.size());
        p("Unaccounted differences but with hotspots same apart from position: " + unaccountedDifferencesButWithHotspotsSameApartFromPosition.size());
//        p("Number of annotations that have same hotspots and use non canonical transcript: " + numberWithSameHotspotsAndUsingNonCanonicalTranscript);
//        p("Number of annotations that have different hotspots and use non canonical transcript: " + numberWithDifferentHotspotsAndUsingNonCanonicalTranscript);
    }

//    @Test
    public void examples()
    {
//        ProteinVariant variant = transval.variationParser().parseExpressionForGene("PLCB4", "M549_G556delinsI");
//        ProteinVariant variant = transval.variationParser().parseExpressionForGene("ROS1", "A1924_I1934del");
        ProteinVariant variant = transval.variationParser().parseVariantForGene("ALK", "D1276_R1279delinsE");
        TransvalVariant tsm = variant.calculateVariant(transval.mRefGenome);
        Assert.assertEquals(6, tsm.hotspots().size());
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
        final StatsForGene statsForGene = new StatsForGene(gene);
        collators.keySet().forEach(annotation ->
        {
            ProteinAnnotationCollator collator = collators.get(annotation);
            VariantStatus comparison = checkVariant(collator);
//            if(collator.mGene.contains("ARID1A"))
//            {
                p( gene + " " + annotation);
//            }
            if(comparison.parsedOk())
            {
                if(comparison.hasProcessingError())
                {
                    p("Processing error for: " + collator);
                }
                else
                {
                    statsForGene.recordResultsForAnnotation(comparison);
                }
            }
            else
            {
                statsForGene.recordNotParsed();
            }
        });
        return statsForGene;
    }

    private VariantStatus checkVariant(ProteinAnnotationCollator collator)
    {
        ProteinVariant variant;
        try
        {
            variant = transval.variationParser().parseVariantForGene(collator.mGene, collator.mAnnotation);
        }
        catch(Exception e)
        {
            return new VariantStatus(collator, e);
        }
        try
        {
            TransvalVariant transvalVariant = variant.calculateVariant(transval.mRefGenome);
            return new VariantStatus(collator, transvalVariant);
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

    private TransvalVariant variant;

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, Exception parseException)
    {
        this.collator = collator;
        mParseException = parseException;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, TransvalVariant variant)
    {
        this.collator = collator;
        this.variant = variant;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator)
    {
        this.collator = collator;
    }

    GeneAnnotation geneAnnotation()
    {
        return new GeneAnnotation(collator.mGene, collator.mAnnotation);
    }

    boolean usesNonCanonicalTranscript()
    {
        return !variant.Transcript.IsCanonical;
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
        // The serve data contain Hotspots such as
        // Hotspot{ref=GCAACATCTCCGAAAGCCAACAAGGAAAT, alt=GG, chromosome=chr7, position=55174785}
        // This has a prefix common to the ref and alt that should be removed
        // and the position should be adjusted accordingly.
        Set<TransvalHotspot> reducedCollatorResults = collator.hotspots.stream().map(this::reduceHotspot).collect(Collectors.toSet());
        // For dels, we only produce one hotspot. If it's in the reduced hotspots, that's good enough.
        if(collator.mAnnotation.endsWith("del"))
        {
            TransvalHotspot computed = variant.hotspots().iterator().next();
//            if(collator.mAnnotation.endsWith("del"))
//            {
//                System.out.println("del");
//            }
            if(reducedCollatorResults.contains(computed))
            {
                return true;
            }
            else
            {
                System.out.println("----------- " + collator.mGene + " " + collator.mAnnotation);
                System.out.println("Computed: " + computed);
                System.out.println("Reduced: " + reducedCollatorResults);
                return false;
            }
        }
        if(collator.mAnnotation.endsWith("dup"))
        {
            return variant.hotspots().containsAll(variant.hotspots());
        }
        final boolean equals = variant.hotspots().containsAll(reducedCollatorResults);
        return equals;
    }

    boolean hotspotsSameModuloPosition()
    {
        Set<FloatingHotspot> variantFloatingHotspots = variant.hotspots().stream().map(FloatingHotspot::new).collect(Collectors.toSet());
        Set<FloatingHotspot> collatorFloatingHotspots = collator.hotspots.stream().map(FloatingHotspot::new).collect(Collectors.toSet());
        return variantFloatingHotspots.equals(collatorFloatingHotspots);
    }

    private TransvalHotspot reduceHotspot(TransvalHotspot hotspot)
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

    public FloatingHotspot(TransvalHotspot hotspot)
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

class StatsForGene
{
    @NotNull
    public final String mGene;

    public int NumberNotParsed = 0;

    public int NumberWithDifferentHotspots = 0;

    public int NumberWithSameHotspots = 0;

    public int NumberWithSameHotspotsThatUseNonCanonicalTranscript = 0;
    public int NumberWithDifferentHotspotsThatUseNonCanonicalTranscript = 0;

    public final List<VariantStatus> AnnotationsWithDifferentHotspots = new ArrayList<>();
    private final List<VariantStatus> AnnotationResults = new ArrayList<>();

    StatsForGene(@NotNull final String gene)
    {
        mGene = gene;
    }

    public void recordResultsForAnnotation(VariantStatus comparison)
    {
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
            if(comparison.usesNonCanonicalTranscript())
            {
                NumberWithDifferentHotspotsThatUseNonCanonicalTranscript++;
            }
        }
        AnnotationResults.add(comparison);
    }

    public List<VariantStatus> getAnnotationsWithDifferentHotspots()
    {
        return AnnotationsWithDifferentHotspots;
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
    Set<TransvalHotspot> hotspots = new HashSet<>();

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
        hotspots.add(new TransvalHotspot(serveItem.Ref, serveItem.Alt, mChromosome, serveItem.Position));
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

