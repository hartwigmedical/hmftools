package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

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
    private final File ensemblDataDir = new File("/Users/timlavers/work/data/v6_0/ref/38/common/ensembl_data");
    private Transval transval;

    @Before
    public void setup() throws IOException
    {
        loadServeItems();
        createTransval();
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
        var map =gson.fromJson(contents, JsonObject.class);
        var records = map.get("records");
        var hotspotsArray = ((JsonObject) records).get("V38").getAsJsonObject().get("knownEvents").getAsJsonObject().get("hotspots").getAsJsonArray();
//        System.out.println(hotspotsArray.size());
//        int numberOfHotspots = hotspotsArray.size();
        for (int i=10_000; i<20_000; i++)
        {
            JsonObject hotspot = hotspotsArray.get(i).getAsJsonObject();
            final String gene = str(hotspot,"gene");
            final String annotation = str(hotspot,"inputProteinAnnotation");
            final String chromosome = str(hotspot,"chromosome");
            final String ref = str(hotspot,"ref");
            final String alt = str(hotspot,"alt");
            final int position = hotspot.get("position").getAsInt();
            if (gene.equals("GLI2") && annotation.equals("R473L"))
            {
                p("GLI2........");
            }
            ServeItem item = new ServeItem(gene, annotation, chromosome, ref, alt, position);
            if (!geneToItems.containsKey(gene))
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

    @Test
    public void check()
    {
        for (String gene : geneToItems.keySet())
        {
            System.out.println("-------------------------");
            System.out.println(gene);
//            System.out.println(geneToItems.get(gene).size());
            checkItemsForGene(gene);
            System.out.println("-------------------------");
        }
    }

    @Test
    public void examples()
    {
        SingleAminoAcidVariant variant = transval.variationParser().parseSingleAminoAcidVariant("GLI2", "R473L");
        TransvalSnvMnv tsm = variant.calculateVariant(transval.mRefGenome);
        Assert.assertEquals(6, tsm.hotspots().size());

    }

    private void checkItemsForGene(String gene)
    {
        Map<String, ProteinAnnotationCollator> collators = new HashMap<>();
        geneToItems.get(gene).forEach(serveItem -> {
            if(!collators.containsKey(serveItem.Annotation))
            {
                collators.put(serveItem.Annotation, new ProteinAnnotationCollator(serveItem));
            }
            collators.get(serveItem.Annotation).addHotspot(serveItem);
        });
        AtomicInteger numberOk = new AtomicInteger(0);
        AtomicInteger numberWithDifferentHotspots = new AtomicInteger(0);
        AtomicInteger numberNotParsed = new AtomicInteger(0);
        collators.keySet().forEach(annotation -> {
            ProteinAnnotationCollator collator = collators.get(annotation);
            VariantStatus comparison = checkVariant(collator);
            if(comparison.parsedOk())
            {
                if(comparison.hasProcessingError())
                {
                    p("Processing error for: " + collator);
                } else
                {
                    if(comparison.hotspotsSame())
                    {
                        //                    System.out.println("Handled OK: " + collator);
                        numberOk.incrementAndGet();
                    }
                    else
                    {
                        p("Hotspots different: " + collator);
                        numberWithDifferentHotspots.incrementAndGet();
                    }
                }
            } else {
//                p("Not parsed: " + collator);
                numberNotParsed.incrementAndGet();
            }
        });
        p("Number not parsed: " + numberNotParsed.get());
        p("Number ok: " + numberOk.get());
        p("Number with different hotspots: " + numberWithDifferentHotspots.get());
    }

    private VariantStatus checkVariant(ProteinAnnotationCollator collator)
    {
        SingleAminoAcidVariant variant;
        try {
            variant = transval.variationParser().parseSingleAminoAcidVariant(collator.mGene, collator.mAnnotation);
        }
        catch(Exception e)
        {
            return new VariantStatus(collator, e);
        }
        try {
            TransvalVariant transvalVariant = variant.calculateVariant(transval.mRefGenome);
            return new VariantStatus(collator, transvalVariant.hotspots());
        } catch(Throwable t){
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
    static VariantStatus withProcessingException(@NotNull final ProteinAnnotationCollator collator, final  @NotNull Throwable e)
    {
        VariantStatus status = new VariantStatus((collator));
        status.mProcessingError = e;
        return status;
    }
    @NotNull
    final ProteinAnnotationCollator collator;

    Exception mParseException = null;

    Throwable mProcessingError = null;

    Set<TransvalHotspot> mHotspotsCalculated = new HashSet<>();

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, Exception parseException)
    {
        this.collator = collator;
        mParseException = parseException;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator, Set<TransvalHotspot> hotspotsCalculated)
    {
        this.collator = collator;
        this.mHotspotsCalculated = hotspotsCalculated;
    }

    VariantStatus(@NotNull final ProteinAnnotationCollator collator)
    {
        this.collator = collator;
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
        return collator.hotspots.equals(mHotspotsCalculated);
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
