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
//                        if (gene.equals("KDM4C") && annotation.equals("G302E"))
//            if (annotation.contains("delins"))
//            {
//                p(annotation);
//            }
            ServeItem item = new ServeItem(gene, annotation, chromosome, ref, alt, position);
            if(!geneToItems.containsKey(gene))
            {
                geneToItems.put(gene, new HashSet<>());
            }
            geneToItems.get(gene).add(item);
        }
    }

    private void loadKnownDifferences()
    {
        knownDifferences = new HashSet<>();
        knownDifferences.add(new GeneAnnotation("KDM4C", "G302E"));
        knownDifferences.add(new GeneAnnotation("KDM4C", "C778F"));
        knownDifferences.add(new GeneAnnotation("KDM4C", "K53M"));
        //        knownDifferences.add(new KnownDifference("", ""));
    }

    private static String str(JsonObject jsonObject, String key)
    {
        return jsonObject.get(key).getAsString();
    }

    @Test
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
        int numberUsingNonCanonicalTranscript = 0;
        SortedSet<GeneAnnotation> unaccountedDifferences = new TreeSet<>();
        for(StatsForGene stats : statsForGenes)
        {
            numberOfGenes++;
            numberNotParsed += stats.NumberNotParsed;
            numberWithSameHotspots += stats.NumberWithSameHotspots;
            numberUsingNonCanonicalTranscript += stats.NumberUsingNonCanonicalTranscript;
            stats.AnnotationsWithDifferentHotspots.forEach(annotation ->
            {
                if(!knownDifferences.contains(annotation))
                {
                    unaccountedDifferences.add(annotation);
                }
            });
        }
        p("UNNACOUNTED DIFFERENCES:");
        p("************************");
//        unaccountedDifferences.forEach(difference ->
//                p(difference.mGene + " " + difference.mAnnotation));
        p("OVERALL STATS");
        p("*************");
        p("Number of genes: " + numberOfGenes);
        p("Number of annotations not parsed: " + numberNotParsed);
        p("Number of annotations with same hotspots: " + numberWithSameHotspots);
        p("Number of annotations with unaccounted differences: " + unaccountedDifferences.size());
        p("Number using non canonical transcript: " + numberUsingNonCanonicalTranscript);
    }

    @Test
    public void examples()
    {
        SingleAminoAcidVariant variant = transval.variationParser().parseSingleAminoAcidVariant("DYRK1A", "D454H");
        TransvalSnvMnv tsm = variant.calculateVariant(transval.mRefGenome);
        Assert.assertEquals(6, tsm.hotspots().size());

        /*

        UTAF1
        SMARCB1
        FGFR1
        FGFR2

        KDM4C
Hotspots different: ProteinAnnotationCollator{mChromosome='chr9', mGene='KDM4C', mAnnotation='G302E'}
Hotspots different: ProteinAnnotationCollator{mChromosome='chr9', mGene='KDM4C', mAnnotation='K53M'}
Hotspots different: ProteinAnnotationCollator{mChromosome='chr9', mGene='KDM4C', mAnnotation='C778F'}
         */
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
            if(comparison.parsedOk())
            {
                if(comparison.hasProcessingError())
                {
                    p("Processing error for: " + collator);
                }
                else
                {
                    if(!comparison.usesCanonicalTranscript())
                    {
                        statsForGene.recordNonCanonicalTranscript();
                    }
                    if(comparison.hotspotsSame())
                    {
                        statsForGene.recordSameHotspots();
                    }
                    else
                    {
                        statsForGene.recordAnnotationWithDifferentHotspots(annotation);
                    }
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
            variant = transval.variationParser().parseSingleAminoAcidVariant(collator.mGene, collator.mAnnotation);
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

    boolean usesCanonicalTranscript()
    {
        return variant.Transcript.IsCanonical;
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
        return variant.hotspots().containsAll(collator.hotspots);
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

    public int NumberUsingNonCanonicalTranscript = 0;

    public final Set<GeneAnnotation> AnnotationsWithDifferentHotspots = new HashSet<>();

    StatsForGene(@NotNull final String gene)
    {
        mGene = gene;
    }

    public void recordNonCanonicalTranscript()
    {
        NumberUsingNonCanonicalTranscript++;
    }

    public void recordNotParsed()
    {
        NumberNotParsed++;
    }

    public void recordSameHotspots()
    {
        NumberWithSameHotspots++;
    }

    public void recordAnnotationWithDifferentHotspots(@NotNull final String annotation)
    {
        AnnotationsWithDifferentHotspots.add(new GeneAnnotation(mGene, annotation));
        NumberWithDifferentHotspots++;
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
